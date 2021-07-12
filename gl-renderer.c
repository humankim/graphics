/*
 * Copyright © 2012 Intel Corporation
 * Copyright © 2015,2019 Collabora, Ltd.
 * Copyright © 2016,2017,2020 NVIDIA Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial
 * portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "config.h"
#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>
#include <linux/input.h>
#include <drm_fourcc.h>
#include <unistd.h>
#include "linux-sync-file.h"
#include "timeline.h"
#include "gl-renderer.h"
#include "gl-renderer-internal.h"
#include "vertex-clipping.h"
#include "linux-dmabuf.h"
#include "linux-dmabuf-unstable-v1-server-protocol.h"
#include "linux-explicit-synchronization.h"
#include "pixel-formats.h"
#include "shared/fd-util.h"
#include "shared/helpers.h"
#include "shared/platform.h"
#include "shared/timespec-util.h"
#include "shared/weston-egl-ext.h"
#include <libweston/lgeFeature.h>
#define GR_GL_VERSION(major, minor) \
	(((uint32_t)(major) << 16) | (uint32_t)(minor))
#define GR_GL_VERSION_INVALID \
	GR_GL_VERSION(0, 0)
#define BUFFER_DAMAGE_COUNT 2
#if HKMC_CCIC
#include <math.h>
#define QT_FAST_BLUR 1
//#define QT_GAUSSIAN_BLUR 0
#define TILT_ENABLE 1
#define REAR_EACH_WIDTH 1280
#if QT_FAST_BLUR
#define FBO_REPAINT_MAX 6
#define DBG_BLUR_LOG 0
#elif QT_GAUSSIAN_BLUR
#define FBO_REPAINT_MAX 2
#define GAUSS_ARRAY_LENGTH 15
#define DBG_BLUR_LOG 1
#endif
#endif // HKMC_CCIC
enum gl_border_status {
	BORDER_STATUS_CLEAN = 0,
	BORDER_TOP_DIRTY = 1 << GL_RENDERER_BORDER_TOP,
	BORDER_LEFT_DIRTY = 1 << GL_RENDERER_BORDER_LEFT,
	BORDER_RIGHT_DIRTY = 1 << GL_RENDERER_BORDER_RIGHT,
	BORDER_BOTTOM_DIRTY = 1 << GL_RENDERER_BORDER_BOTTOM,
	BORDER_ALL_DIRTY = 0xf,
	BORDER_SIZE_CHANGED = 0x10
};
struct gl_border_image {
	GLuint tex;
	int32_t width, height;
	int32_t tex_width;
	void *data;
};
struct gl_output_state {
	EGLSurface egl_surface;
	pixman_region32_t buffer_damage[BUFFER_DAMAGE_COUNT];
	int buffer_damage_index;
	enum gl_border_status border_damage[BUFFER_DAMAGE_COUNT];
	struct gl_border_image borders[4];
	enum gl_border_status border_status;
#if 1 // nvidia
	EGLStreamKHR egl_stream;
#endif
	struct weston_matrix output_matrix;
	EGLSyncKHR begin_render_sync, end_render_sync;
	/* struct timeline_render_point::link */
	struct wl_list timeline_render_point_list;
#if 1 // nvidia
	float gamma[3];  /* For exponential modulation over normalized RGB values */
	float shade[3];  /* Direct fraction multiplication to RGB value */
#endif
#if HKMC_CCIC
#if defined(QT_FAST_BLUR) || defined(QT_GAUSSIAN_BLUR)
	GLuint blur_fbo[FBO_REPAINT_MAX];
	GLuint blur_tex[FBO_REPAINT_MAX];
#endif
#endif
};
enum buffer_type {
	BUFFER_TYPE_NULL,
	BUFFER_TYPE_SOLID, /* internal solid color surfaces without a buffer */
	BUFFER_TYPE_SHM,
	BUFFER_TYPE_EGL
};
struct gl_renderer;
struct egl_image {
	struct gl_renderer *renderer;
	EGLImageKHR image;
	int refcount;
};
enum import_type {
	IMPORT_TYPE_INVALID,
	IMPORT_TYPE_DIRECT,
	IMPORT_TYPE_GL_CONVERSION
};
struct dmabuf_image {
	struct linux_dmabuf_buffer *dmabuf;
	int num_images;
	struct egl_image *images[3];
	struct wl_list link;
	enum import_type import_type;
	GLenum target;
	struct gl_shader *shader;
};
struct yuv_plane_descriptor {
	int width_divisor;
	int height_divisor;
	uint32_t format;
	int plane_index;
};
enum texture_type {
	TEXTURE_Y_XUXV_WL,
	TEXTURE_Y_UV_WL,
	TEXTURE_Y_U_V_WL,
	TEXTURE_XYUV_WL
};
struct yuv_format_descriptor {
	uint32_t format;
	int input_planes;
	int output_planes;
	enum texture_type texture_type;
	struct yuv_plane_descriptor plane[4];
};
struct gl_surface_state {
	GLfloat color[4];
	struct gl_shader *shader;
	GLuint textures[3];
	int num_textures;
	bool needs_full_upload;
	pixman_region32_t texture_damage;
	/* These are only used by SHM surfaces to detect when we need
	 * to do a full upload to specify a new internal texture
	 * format */
	GLenum gl_format[3];
	GLenum gl_pixel_type;
	struct egl_image* images[3];
	GLenum target;
	int num_images;
	struct weston_buffer_reference buffer_ref;
	struct weston_buffer_release_reference buffer_release_ref;
	enum buffer_type buffer_type;
	int pitch; /* in pixels */
	int height; /* in pixels */
	bool y_inverted;
	bool direct_display;
	/* Extension needed for SHM YUV texture */
	int offset[3]; /* offset per plane */
	int hsub[3];  /* horizontal subsampling per plane */
	int vsub[3];  /* vertical subsampling per plane */
#if 1 // nvidia
	EGLStreamKHR egl_stream;
	bool new_stream;
#endif
	struct weston_surface *surface;
	/* Whether this surface was used in the current output repaint.
	   Used only in the context of a gl_renderer_repaint_output call. */
	bool used_in_output_repaint;
	struct wl_listener surface_destroy_listener;
	struct wl_listener renderer_destroy_listener;
};
enum timeline_render_point_type {
	TIMELINE_RENDER_POINT_TYPE_BEGIN,
	TIMELINE_RENDER_POINT_TYPE_END
};
struct timeline_render_point {
	struct wl_list link; /* gl_output_state::timeline_render_point_list */
	enum timeline_render_point_type type;
	int fd;
	struct weston_output *output;
	struct wl_event_source *event_source;
};
#if HKMC_CCIC
#if defined(TILT_ENABLE) || defined(QT_FAST_BLUR) || defined(QT_GAUSSIAN_BLUR)
bool is_rear_screen(char *output_name);
#endif
#if defined(QT_FAST_BLUR) || defined(QT_GAUSSIAN_BLUR)
void convert_tcoord(struct weston_output *output, pixman_box32_t *region, GLfloat *tcoord);
void scaling_blur_coord(struct weston_output *output, pixman_box32_t *rects, GLfloat *verts);
void create_fbo(struct weston_output *output, struct gl_output_state *go);
void bind_fbo(struct weston_view *ev, struct weston_output *output);
#endif
#if QT_FAST_BLUR
GLfloat get_weight(GLfloat v);
void repaint_blur(struct weston_view *ev, struct weston_output *output, pixman_region32_t *damage);
#elif QT_GAUSSIAN_BLUR
struct gauss_sample
{
	char* name;
	double pos;
	double weight;
};
char *itoa(int num, char *str);
void set_gauss_sample(struct gauss_sample *sample, char* n, double p, double w);
void set_gauss_properties(struct weston_view *ev, struct weston_output *output, struct gl_shader *shader,
			    bool is_horizontal, float *spread, uint32_t *samples, float (* dirstep)[2]);
void repaint_gaussian_blur(struct weston_view *ev, struct weston_output *output,
				pixman_region32_t *damage, bool is_horizontal);
#endif
#endif // HKMC_CCIC
static PFNEGLGETPLATFORMDISPLAYEXTPROC get_platform_display = NULL;
static inline const char *
dump_format(uint32_t format, char out[4])
{
#if BYTE_ORDER == BIG_ENDIAN
	format = __builtin_bswap32(format);
#endif
	memcpy(out, &format, 4);
	return out;
}
static inline struct gl_output_state *
get_output_state(struct weston_output *output)
{
	return (struct gl_output_state *)output->renderer_state;
}
static int
gl_renderer_create_surface(struct weston_surface *surface);
static inline struct gl_surface_state *
get_surface_state(struct weston_surface *surface)
{
	if (!surface->renderer_state)
		gl_renderer_create_surface(surface);
	return (struct gl_surface_state *)surface->renderer_state;
}
static void
timeline_render_point_destroy(struct timeline_render_point *trp)
{
	wl_list_remove(&trp->link);
	wl_event_source_remove(trp->event_source);
	close(trp->fd);
	free(trp);
}
static int
timeline_render_point_handler(int fd, uint32_t mask, void *data)
{
	struct timeline_render_point *trp = data;
	const char *tp_name = trp->type == TIMELINE_RENDER_POINT_TYPE_BEGIN ?
			      "renderer_gpu_begin" : "renderer_gpu_end";
	if (mask & WL_EVENT_READABLE) {
		struct timespec tspec = { 0 };
		if (weston_linux_sync_file_read_timestamp(trp->fd,
							  &tspec) == 0) {
			TL_POINT(trp->output->compositor, tp_name, TLP_GPU(&tspec),
				 TLP_OUTPUT(trp->output), TLP_END);
		}
	}
	timeline_render_point_destroy(trp);
	return 0;
}
static EGLSyncKHR
create_render_sync(struct gl_renderer *gr)
{
	static const EGLint attribs[] = { EGL_NONE };
	if (!gr->has_native_fence_sync)
		return EGL_NO_SYNC_KHR;
	return gr->create_sync(gr->egl_display, EGL_SYNC_NATIVE_FENCE_ANDROID,
			       attribs);
}
static void
timeline_submit_render_sync(struct gl_renderer *gr,
			    struct weston_compositor *ec,
			    struct weston_output *output,
			    EGLSyncKHR sync,
			    enum timeline_render_point_type type)
{
	struct gl_output_state *go;
	struct wl_event_loop *loop;
	int fd;
	struct timeline_render_point *trp;
	if (!weston_log_scope_is_enabled(ec->timeline) ||
	    !gr->has_native_fence_sync ||
	    sync == EGL_NO_SYNC_KHR)
		return;
	go = get_output_state(output);
	loop = wl_display_get_event_loop(ec->wl_display);
	fd = gr->dup_native_fence_fd(gr->egl_display, sync);
	if (fd == EGL_NO_NATIVE_FENCE_FD_ANDROID)
		return;
	trp = zalloc(sizeof *trp);
	if (trp == NULL) {
		close(fd);
		return;
	}
	trp->type = type;
	trp->fd = fd;
	trp->output = output;
	trp->event_source = wl_event_loop_add_fd(loop, fd,
						 WL_EVENT_READABLE,
						 timeline_render_point_handler,
						 trp);
	wl_list_insert(&go->timeline_render_point_list, &trp->link);
}
static struct egl_image*
egl_image_create(struct gl_renderer *gr, EGLenum target,
		 EGLClientBuffer buffer, const EGLint *attribs)
{
	struct egl_image *img;
	img = zalloc(sizeof *img);
	img->renderer = gr;
	img->refcount = 1;
	img->image = gr->create_image(gr->egl_display, EGL_NO_CONTEXT,
				      target, buffer, attribs);
	if (img->image == EGL_NO_IMAGE_KHR) {
		free(img);
		return NULL;
	}
	return img;
}
static struct egl_image*
egl_image_ref(struct egl_image *image)
{
	image->refcount++;
	return image;
}
static int
egl_image_unref(struct egl_image *image)
{
	struct gl_renderer *gr = image->renderer;
	assert(image->refcount > 0);
	image->refcount--;
	if (image->refcount > 0)
		return image->refcount;
	gr->destroy_image(gr->egl_display, image->image);
	free(image);
	return 0;
}
static struct dmabuf_image*
dmabuf_image_create(void)
{
	struct dmabuf_image *img;
	img = zalloc(sizeof *img);
	wl_list_init(&img->link);
	return img;
}
static void
dmabuf_image_destroy(struct dmabuf_image *image)
{
	int i;
	for (i = 0; i < image->num_images; ++i)
		egl_image_unref(image->images[i]);
	if (image->dmabuf)
		linux_dmabuf_buffer_set_user_data(image->dmabuf, NULL, NULL);
	wl_list_remove(&image->link);
	free(image);
}
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) > (b)) ? (b) : (a))
/*
 * Compute the boundary vertices of the intersection of the global coordinate
 * aligned rectangle 'rect', and an arbitrary quadrilateral produced from
 * 'surf_rect' when transformed from surface coordinates into global coordinates.
 * The vertices are written to 'ex' and 'ey', and the return value is the
 * number of vertices. Vertices are produced in clockwise winding order.
 * Guarantees to produce either zero vertices, or 3-8 vertices with non-zero
 * polygon area.
 */
static int
calculate_edges(struct weston_view *ev, pixman_box32_t *rect,
		pixman_box32_t *surf_rect, GLfloat *ex, GLfloat *ey)
{
	struct clip_context ctx;
	int i, n;
	GLfloat min_x, max_x, min_y, max_y;
	struct polygon8 surf = {
		{ surf_rect->x1, surf_rect->x2, surf_rect->x2, surf_rect->x1 },
		{ surf_rect->y1, surf_rect->y1, surf_rect->y2, surf_rect->y2 },
		4
	};
	ctx.clip.x1 = rect->x1;
	ctx.clip.y1 = rect->y1;
	ctx.clip.x2 = rect->x2;
	ctx.clip.y2 = rect->y2;
	/* transform surface to screen space: */
	for (i = 0; i < surf.n; i++)
		weston_view_to_global_float(ev, surf.x[i], surf.y[i],
					    &surf.x[i], &surf.y[i]);
	/* find bounding box: */
	min_x = max_x = surf.x[0];
	min_y = max_y = surf.y[0];
	for (i = 1; i < surf.n; i++) {
		min_x = min(min_x, surf.x[i]);
		max_x = max(max_x, surf.x[i]);
		min_y = min(min_y, surf.y[i]);
		max_y = max(max_y, surf.y[i]);
	}
	/* First, simple bounding box check to discard early transformed
	 * surface rects that do not intersect with the clip region:
	 */
	if ((min_x >= ctx.clip.x2) || (max_x <= ctx.clip.x1) ||
	    (min_y >= ctx.clip.y2) || (max_y <= ctx.clip.y1))
		return 0;
	/* Simple case, bounding box edges are parallel to surface edges,
	 * there will be only four edges.  We just need to clip the surface
	 * vertices to the clip rect bounds:
	 */
	if (!ev->transform.enabled)
		return clip_simple(&ctx, &surf, ex, ey);
	/* Transformed case: use a general polygon clipping algorithm to
	 * clip the surface rectangle with each side of 'rect'.
	 * The algorithm is Sutherland-Hodgman, as explained in
	 * http://www.codeguru.com/cpp/misc/misc/graphics/article.php/c8965/Polygon-Clipping.htm
	 * but without looking at any of that code.
	 */
	n = clip_transformed(&ctx, &surf, ex, ey);
	if (n < 3)
		return 0;
	return n;
}
static bool
merge_down(pixman_box32_t *a, pixman_box32_t *b, pixman_box32_t *merge)
{
	if (a->x1 == b->x1 && a->x2 == b->x2 && a->y1 == b->y2) {
		merge->x1 = a->x1;
		merge->x2 = a->x2;
		merge->y1 = b->y1;
		merge->y2 = a->y2;
		return true;
	}
	return false;
}
static int
compress_bands(pixman_box32_t *inrects, int nrects,
		   pixman_box32_t **outrects)
{
	bool merged = false;
	pixman_box32_t *out, merge_rect;
	int i, j, nout;
	if (!nrects) {
		*outrects = NULL;
		return 0;
	}
	/* nrects is an upper bound - we're not too worried about
	 * allocating a little extra
	 */
	out = malloc(sizeof(pixman_box32_t) * nrects);
	out[0] = inrects[0];
	nout = 1;
	for (i = 1; i < nrects; i++) {
		for (j = 0; j < nout; j++) {
			merged = merge_down(&inrects[i], &out[j], &merge_rect);
			if (merged) {
				out[j] = merge_rect;
				break;
			}
		}
		if (!merged) {
			out[nout] = inrects[i];
			nout++;
		}
	}
	*outrects = out;
	return nout;
}
static int
texture_region(struct weston_view *ev, pixman_region32_t *region,
		pixman_region32_t *surf_region)
{
	struct gl_surface_state *gs = get_surface_state(ev->surface);
	struct weston_compositor *ec = ev->surface->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	GLfloat *v, inv_width, inv_height;
	unsigned int *vtxcnt, nvtx = 0;
	pixman_box32_t *rects, *surf_rects;
	pixman_box32_t *raw_rects;
	int i, j, k, nrects, nsurf, raw_nrects;
	bool used_band_compression;
	raw_rects = pixman_region32_rectangles(region, &raw_nrects);
	surf_rects = pixman_region32_rectangles(surf_region, &nsurf);
	if (raw_nrects < 4) {
		used_band_compression = false;
		nrects = raw_nrects;
		rects = raw_rects;
	} else {
		nrects = compress_bands(raw_rects, raw_nrects, &rects);
		used_band_compression = true;
	}
	/* worst case we can have 8 vertices per rect (ie. clipped into
	 * an octagon):
	 */
	v = wl_array_add(&gr->vertices, nrects * nsurf * 8 * 4 * sizeof *v);
	vtxcnt = wl_array_add(&gr->vtxcnt, nrects * nsurf * sizeof *vtxcnt);
	inv_width = 1.0 / gs->pitch;
        inv_height = 1.0 / gs->height;
	for (i = 0; i < nrects; i++) {
		pixman_box32_t *rect = &rects[i];
		for (j = 0; j < nsurf; j++) {
			pixman_box32_t *surf_rect = &surf_rects[j];
			GLfloat sx, sy, bx, by;
			GLfloat ex[8], ey[8];          /* edge points in screen space */
			int n;
			/* The transformed surface, after clipping to the clip region,
			 * can have as many as eight sides, emitted as a triangle-fan.
			 * The first vertex in the triangle fan can be chosen arbitrarily,
			 * since the area is guaranteed to be convex.
			 *
			 * If a corner of the transformed surface falls outside of the
			 * clip region, instead of emitting one vertex for the corner
			 * of the surface, up to two are emitted for two corresponding
			 * intersection point(s) between the surface and the clip region.
			 *
			 * To do this, we first calculate the (up to eight) points that
			 * form the intersection of the clip rect and the transformed
			 * surface.
			 */
			n = calculate_edges(ev, rect, surf_rect, ex, ey);
			if (n < 3)
				continue;
			/* emit edge points: */
			for (k = 0; k < n; k++) {
				weston_view_from_global_float(ev, ex[k], ey[k],
							      &sx, &sy);
				/* position: */
				*(v++) = ex[k];
				*(v++) = ey[k];
				/* texcoord: */
				weston_surface_to_buffer_float(ev->surface,
							       sx, sy,
							       &bx, &by);
				*(v++) = bx * inv_width;
				if (gs->y_inverted) {
					*(v++) = by * inv_height;
				} else {
					*(v++) = (gs->height - by) * inv_height;
				}
			}
			vtxcnt[nvtx++] = n;
		}
	}
	if (used_band_compression)
		free(rects);
	return nvtx;
}
static void
triangle_fan_debug(struct weston_view *view, int first, int count)
{
	struct weston_compositor *compositor = view->surface->compositor;
	struct gl_renderer *gr = get_renderer(compositor);
	int i;
	GLushort *buffer;
	GLushort *index;
	int nelems;
	static int color_idx = 0;
	static const GLfloat color[][4] = {
			{ 1.0, 0.0, 0.0, 1.0 },
			{ 0.0, 1.0, 0.0, 1.0 },
			{ 0.0, 0.0, 1.0, 1.0 },
			{ 1.0, 1.0, 1.0, 1.0 },
	};
	nelems = (count - 1 + count - 2) * 2;
	buffer = malloc(sizeof(GLushort) * nelems);
	index = buffer;
	for (i = 1; i < count; i++) {
		*index++ = first;
		*index++ = first + i;
	}
	for (i = 2; i < count; i++) {
		*index++ = first + i - 1;
		*index++ = first + i;
	}
	glUseProgram(gr->solid_shader.program);
	glUniform4fv(gr->solid_shader.color_uniform, 1,
			color[color_idx++ % ARRAY_LENGTH(color)]);
	glDrawElements(GL_LINES, nelems, GL_UNSIGNED_SHORT, buffer);
	glUseProgram(gr->current_shader->program);
	free(buffer);
}
static void
repaint_region(struct weston_view *ev, pixman_region32_t *region,
		pixman_region32_t *surf_region)
{
	struct weston_compositor *ec = ev->surface->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	GLfloat *v;
	unsigned int *vtxcnt;
	int i, first, nfans;
	/* The final region to be painted is the intersection of
	 * 'region' and 'surf_region'. However, 'region' is in the global
	 * coordinates, and 'surf_region' is in the surface-local
	 * coordinates. texture_region() will iterate over all pairs of
	 * rectangles from both regions, compute the intersection
	 * polygon for each pair, and store it as a triangle fan if
	 * it has a non-zero area (at least 3 vertices, actually).
	 */
	nfans = texture_region(ev, region, surf_region);
	v = gr->vertices.data;
	vtxcnt = gr->vtxcnt.data;
	/* position: */
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof *v, &v[0]);
	glEnableVertexAttribArray(0);
	/* texcoord: */
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof *v, &v[2]);
	glEnableVertexAttribArray(1);
	for (i = 0, first = 0; i < nfans; i++) {
		glDrawArrays(GL_TRIANGLE_FAN, first, vtxcnt[i]);
		if (gr->fan_debug)
			triangle_fan_debug(ev, first, vtxcnt[i]);
		first += vtxcnt[i];
	}
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	gr->vertices.size = 0;
	gr->vtxcnt.size = 0;
}
static int
use_output(struct weston_output *output)
{
	static int errored;
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	EGLBoolean ret;
	ret = eglMakeCurrent(gr->egl_display, go->egl_surface,
			     go->egl_surface, gr->egl_context);
	if (ret == EGL_FALSE) {
		if (errored)
			return -1;
		errored = 1;
		weston_log("Failed to make EGL context current.\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	return 0;
}
static int
shader_init(struct gl_shader *shader, struct gl_renderer *gr,
		   const char *vertex_source, const char *fragment_source);
static void
use_shader(struct gl_renderer *gr, struct gl_shader *shader)
{
	if (!shader->program) {
		int ret;
		ret =  shader_init(shader, gr,
				   shader->vertex_source,
				   shader->fragment_source);
		if (ret < 0)
			weston_log("warning: failed to compile shader\n");
	}
	if (gr->current_shader == shader)
		return;
	glUseProgram(shader->program);
	gr->current_shader = shader;
}
#if HKMC_CCIC
#if TILT_ENABLE
static void init_output_matrix(struct weston_output *output, bool tilt_on, int32_t x_offset) {
	struct gl_output_state *go = get_output_state(output);
	go->output_matrix = output->matrix;
	int scaling_factor = 1;
	if (tilt_on) {
		scaling_factor = 2;
	}
	glViewport(
		x_offset + go->borders[GL_RENDERER_BORDER_LEFT].width,
		go->borders[GL_RENDERER_BORDER_BOTTOM].height,
		output->current_mode->width / scaling_factor,
		output->current_mode->height);
	weston_matrix_translate(&go->output_matrix,
		-(x_offset + (output->current_mode->width) / (2.0 * scaling_factor)),
		-(output->current_mode->height / 2.0), 0);
	weston_matrix_scale(&go->output_matrix,
		(2.0 * scaling_factor) / output->current_mode->width,
		-2.0 / output->current_mode->height, 1);
}
#endif
#if QT_GAUSSIAN_BLUR
char *itoa(int num, char *str)
{
	if(str == NULL)
		return NULL;
	sprintf(str, "%d", num);
	return str;
}
static double gaussian(double x, double d)
{
	return exp(- x * x / (2 * d * d));
}
void set_gauss_sample(struct gauss_sample *sample, char* n, double p, double w)
{
	sample->name = n;
	sample->pos = p;
	sample->weight = w;
}
static void build_gauss_sample_points(struct gl_shader *shader, int samples, uint32_t radius, double deviation)
{
	char s_index[3];
	uint32_t i=0;
	double p0, p1, w0, w1, w, wl, wc;
	double sample_pos;
	char s[10];
	struct gauss_sample *p = malloc(samples * sizeof(struct gauss_sample));
	if ((samples % 2) == 1) {
		set_gauss_sample(&p[radius], "tC", 0.0, 1.0);
		glUniform2f(shader->tC_uniform, 0.0, 1.0);
		for (i=0; i<radius; ++i) {
			p0 = (i + 1) * 2 - 1;
			p1 = (i + 1) * 2;
			w0 = gaussian(p0, deviation);
			w1 = gaussian(p1, deviation);
			w = w0 + w1;
			sample_pos = (p0 * w0 + p1 * w1) / w;
			if (sample_pos <= 0) {
				sample_pos = 0;
				w = 0;
			}
			itoa(i, s_index);
			if(i==0) strcpy(s_index, "0");
			strcpy(s, "tL");
			set_gauss_sample(&p[radius - i - 1], strcat(s,s_index), sample_pos, w);
			glUniform2f(shader->tL_uniform[i], sample_pos, w);
			strcpy(s, "tR");
			set_gauss_sample(&p[radius + i + 1], strcat(s,s_index), -sample_pos, w);
			glUniform2f(shader->tR_uniform[i], -sample_pos, w);
		}
	} else {
		{
		wl = gaussian(-1.0, deviation);
		wc = gaussian(0.0, deviation);
		w = wl + wc;
		set_gauss_sample(&p[radius], "tL0", -1.0 * wl / w, w);
		glUniform2f(shader->tL_uniform[0], -1.0 * wl / w, w);
		set_gauss_sample(&p[radius+1], "tR0", 1.0, wl);
		glUniform2f(shader->tR_uniform[0], 1.0, wl);
		}
		for (i=0; i<radius; ++i) {
			p0 = (i + 1) * 2;
			p1 = (i + 1) * 2 + 1;
			w0 = gaussian(p0, deviation);
			w1 = gaussian(p1, deviation);
			w = w0 + w1;
			sample_pos = (p0 * w0 + p1 * w1) / w;
			if (sample_pos <= 0) {
				sample_pos = 0;
				w = 0;
			}
			itoa(i + 1, s_index);
			strcpy(s, "tL");
			set_gauss_sample(&p[radius - i - 1], strcat(s, s_index), sample_pos, w);
			glUniform2f(shader->tL_uniform[i+1], sample_pos, w);
			strcpy(s, "tR");
			set_gauss_sample(&p[radius + i + 2], strcat(s, s_index), -sample_pos, w);
			glUniform2f(shader->tR_uniform[i+1], -sample_pos, w);
		}
	}
	free(p);
}
#endif
#endif
static void
shader_uniforms(struct gl_shader *shader,
		struct weston_view *view,
		struct weston_output *output)
{
	int i;
	struct gl_surface_state *gs = get_surface_state(view->surface);
	struct gl_output_state *go = get_output_state(output);
#if HKMC_CCIC
	struct weston_matrix tilt_matrix;
	weston_matrix_init(&tilt_matrix);
#if TILT_ENABLE
	if (view->tilt_info.visibility) {
		float rotation_origin_x = (float)wl_fixed_to_double(view->tilt_info.tilt_rotation_origin_x);
		float rotation_origin_y = (float)wl_fixed_to_double(view->tilt_info.tilt_rotation_origin_y);
		float rotation_angle = (float)wl_fixed_to_double(view->tilt_info.tilt_rotation_angle);
		float scale_x = (float)wl_fixed_to_double(view->tilt_info.tilt_scale_x);
		float scale_y = (float)wl_fixed_to_double(view->tilt_info.tilt_scale_y);
		float move_y = (float)wl_fixed_to_double(view->tilt_info.tilt_move_y);
		int32_t x_offset = (int32_t)wl_fixed_to_double(view->tilt_info.tilt_offset_x);
		if (is_rear_screen(output->name)) {
			init_output_matrix(output, true, x_offset);
		}
		if (rotation_angle >= 0) {
			float rotation_origin_mat_x = (2.0 / output->current_mode->width) * rotation_origin_x;
			float rotation_origin_mat_y = (2.0 / output->current_mode->height) * rotation_origin_y;
			float move_mat_y = (2.0 / output->current_mode->height) * move_y;
			weston_matrix_scale(&tilt_matrix, scale_x, scale_y, 1);
			weston_matrix_translate(&tilt_matrix, -rotation_origin_mat_x, -rotation_origin_mat_y, 0);
			weston_matrix_rotate_x(&tilt_matrix, rotation_angle);
			weston_matrix_translate(&tilt_matrix, rotation_origin_mat_x, rotation_origin_mat_y, 0);
			weston_matrix_translate(&tilt_matrix, 0, move_mat_y, 0);
			weston_matrix_frustum(&tilt_matrix, 0.0, 1.0);
		}
	} else if (is_rear_screen(output->name)) {
		init_output_matrix(output, false, 0);
	}
#endif
#endif // HKMC_CCIC
#if HKMC_CCIC
	glUniformMatrix4fv(shader->tilt_uniform, 1, GL_FALSE, tilt_matrix.d);
	glUniformMatrix4fv(shader->proj_uniform, 1, GL_FALSE, go->output_matrix.d);
#else
	glUniformMatrix4fv(shader->proj_uniform,
			   1, GL_FALSE, go->output_matrix.d);
#endif
	glUniformMatrix4fv(shader->proj_uniform,
			   1, GL_FALSE, go->output_matrix.d);
	glUniform4fv(shader->color_uniform, 1, gs->color);
	glUniform1f(shader->alpha_uniform, view->alpha);
#if 1 // nvidia
	glUniform3f(shader->gamma_uniform, go->gamma[0], go->gamma[1], go->gamma[2]);
	glUniform3f(shader->shade_uniform, go->shade[0], go->shade[1], go->shade[2]);
#endif
	for (i = 0; i < gs->num_textures; i++)
		glUniform1i(shader->tex_uniforms[i], i);
#if HKMC_CCIC
	glUniform2i(shader->output_size_uniform, output->current_mode->width, output->current_mode->height);
#endif
}
#if HKMC_CCIC
static bool
has_mask_region(struct weston_view *ev) {
	if (ev->mask_region.x1 > 0 || ev->mask_region.x2 > 0 ||
		ev->mask_region.y1 > 0 || ev->mask_region.y2 > 0) {
		return true;
	}
	return false;
}
#endif
static int
ensure_surface_buffer_is_ready(struct gl_renderer *gr,
			       struct gl_surface_state *gs)
{
	EGLint attribs[] = {
		EGL_SYNC_NATIVE_FENCE_FD_ANDROID,
		-1,
		EGL_NONE
	};
	struct weston_surface *surface = gs->surface;
	struct weston_buffer *buffer = gs->buffer_ref.buffer;
	EGLSyncKHR sync;
	EGLint wait_ret;
	EGLint destroy_ret;
	if (!buffer)
		return 0;
	if (surface->acquire_fence_fd < 0)
		return 0;
	/* We should only get a fence if we support EGLSyncKHR, since
	 * we don't advertise the explicit sync protocol otherwise. */
	assert(gr->has_native_fence_sync);
	/* We should only get a fence for non-SHM buffers, since surface
	 * commit would have failed otherwise. */
	assert(wl_shm_buffer_get(buffer->resource) == NULL);
	attribs[1] = dup(surface->acquire_fence_fd);
	if (attribs[1] == -1) {
		linux_explicit_synchronization_send_server_error(
			gs->surface->synchronization_resource,
			"Failed to dup acquire fence");
		return -1;
	}
	sync = gr->create_sync(gr->egl_display,
			       EGL_SYNC_NATIVE_FENCE_ANDROID,
			       attribs);
	if (sync == EGL_NO_SYNC_KHR) {
		linux_explicit_synchronization_send_server_error(
			gs->surface->synchronization_resource,
			"Failed to create EGLSyncKHR object");
		close(attribs[1]);
		return -1;
	}
	wait_ret = gr->wait_sync(gr->egl_display, sync, 0);
	if (wait_ret == EGL_FALSE) {
		linux_explicit_synchronization_send_server_error(
			gs->surface->synchronization_resource,
			"Failed to wait on EGLSyncKHR object");
		/* Continue to try to destroy the sync object. */
	}
	destroy_ret = gr->destroy_sync(gr->egl_display, sync);
	if (destroy_ret == EGL_FALSE) {
		linux_explicit_synchronization_send_server_error(
			gs->surface->synchronization_resource,
			"Failed to destroy on EGLSyncKHR object");
	}
	return (wait_ret == EGL_TRUE && destroy_ret == EGL_TRUE) ? 0 : -1;
}
 /* Checks if a view needs to be censored on an output
  * Checks for 2 types of censor requirements
  * - recording_censor: Censor protected view when a
  *   protected view is captured.
  * - unprotected_censor: Censor regions of protected views
  *   when displayed on an output which has lower protection capability.
  * Returns the originally stored gl_shader if content censoring is required,
  * NULL otherwise.
  */
static struct gl_shader *
setup_censor_overrides(struct weston_output *output,
		       struct weston_view *ev)
{
	struct gl_shader *replaced_shader = NULL;
	struct weston_compositor *ec = ev->surface->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(ev->surface);
	bool recording_censor =
		(output->disable_planes > 0) &&
		(ev->surface->desired_protection > WESTON_HDCP_DISABLE);
	bool unprotected_censor =
		(ev->surface->desired_protection > output->current_protection);
	if (gs->direct_display) {
		gs->color[0] = 0.40;
		gs->color[1] = 0.0;
		gs->color[2] = 0.0;
		gs->color[3] = 1.0;
		gs->shader = &gr->solid_shader;
		return gs->shader;
	}
	/* When not in enforced mode, the client is notified of the protection */
	/* change, so content censoring is not required */
	if (ev->surface->protection_mode !=
	    WESTON_SURFACE_PROTECTION_MODE_ENFORCED)
		return NULL;
	if (recording_censor || unprotected_censor) {
		replaced_shader = gs->shader;
		gs->color[0] = 0.40;
		gs->color[1] = 0.0;
		gs->color[2] = 0.0;
		gs->color[3] = 1.0;
		gs->shader = &gr->solid_shader;
	}
	return replaced_shader;
}
static void
draw_view(struct weston_view *ev, struct weston_output *output,
	  pixman_region32_t *damage) /* in global coordinates */
{
	struct weston_compositor *ec = ev->surface->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(ev->surface);
	/* repaint bounding region in global coordinates: */
	pixman_region32_t repaint;
	/* opaque region in surface coordinates: */
	pixman_region32_t surface_opaque;
	/* non-opaque region in surface coordinates: */
	pixman_region32_t surface_blend;
	GLint filter;
#if HKMC_CCIC
	/* mask region in surface coordinates: */
	pixman_region32_t mask_region;
#endif
	int i;
	struct gl_shader *replaced_shader = NULL;
	/* In case of a runtime switch of renderers, we may not have received
	 * an attach for this surface since the switch. In that case we don't
	 * have a valid buffer or a proper shader set up so skip rendering. */
	if (!gs->shader && !gs->direct_display)
		return;
	pixman_region32_init(&repaint);
	pixman_region32_intersect(&repaint,
				  &ev->transform.boundingbox, damage);
	pixman_region32_subtract(&repaint, &repaint, &ev->clip);
	if (!pixman_region32_not_empty(&repaint))
		goto out;
	if (ensure_surface_buffer_is_ready(gr, gs) < 0)
		goto out;
	replaced_shader = setup_censor_overrides(output, ev);
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	if (gr->fan_debug) {
		use_shader(gr, &gr->solid_shader);
		shader_uniforms(&gr->solid_shader, ev, output);
	}
	use_shader(gr, gs->shader);
	shader_uniforms(gs->shader, ev, output);
#if HKMC_CCIC
	if (!output->enable_fbo) {
		if (strcmp(output->name, "HDMI-A-1")==0 ||
                strcmp(output->name, "HDMI-A-2")==0) {
			output->enable_fbo = true;
			create_fbo(output, get_output_state(output));
		}
	}
#if QT_FAST_BLUR
	if (output->enable_fbo) {
		bind_fbo(ev, output);
	}
#elif QT_GAUSSIAN_BLUR
	if (output->enable_fbo) {
		bind_fbo(ev, output);
	}
#endif
#endif
	if (ev->transform.enabled || output->zoom.active ||
	    output->current_scale != ev->surface->buffer_viewport.buffer.scale)
		filter = GL_LINEAR;
	else
		filter = GL_NEAREST;
	for (i = 0; i < gs->num_textures; i++) {
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(gs->target, gs->textures[i]);
		glTexParameteri(gs->target, GL_TEXTURE_MIN_FILTER, filter);
		glTexParameteri(gs->target, GL_TEXTURE_MAG_FILTER, filter);
	}
	/* blended region is whole surface minus opaque region: */
	pixman_region32_init_rect(&surface_blend, 0, 0,
				  ev->surface->width, ev->surface->height);
	if (ev->geometry.scissor_enabled)
		pixman_region32_intersect(&surface_blend, &surface_blend,
					  &ev->geometry.scissor);
	pixman_region32_subtract(&surface_blend, &surface_blend,
				 &ev->surface->opaque);
	/* XXX: Should we be using ev->transform.opaque here? */
	pixman_region32_init(&surface_opaque);
	if (ev->geometry.scissor_enabled)
		pixman_region32_intersect(&surface_opaque,
					  &ev->surface->opaque,
					  &ev->geometry.scissor);
	else
		pixman_region32_copy(&surface_opaque, &ev->surface->opaque);
	if (pixman_region32_not_empty(&surface_opaque)) {
		if (gs->shader == &gr->texture_shader_rgba) {
			/* Special case for RGBA textures with possibly
			 * bad data in alpha channel: use the shader
			 * that forces texture alpha = 1.0.
			 * Xwayland surfaces need this.
			 */
			use_shader(gr, &gr->texture_shader_rgbx);
			shader_uniforms(&gr->texture_shader_rgbx, ev, output);
		}
		if (ev->alpha < 1.0)
			glEnable(GL_BLEND);
		else
			glDisable(GL_BLEND);
		repaint_region(ev, &repaint, &surface_opaque);
		gs->used_in_output_repaint = true;
	}
#if HKMC_CCIC
	if (has_mask_region(ev)) {
		pixman_region32_init_rect(&mask_region, ev->mask_region.x1 + output->x,
					  ev->mask_region.y1 + output->y, ev->mask_region.x2,
					  ev->mask_region.y2);
		pixman_region32_subtract(&repaint, &repaint, &mask_region);
	}
#endif
	if (pixman_region32_not_empty(&surface_blend)) {
		use_shader(gr, gs->shader);
		glEnable(GL_BLEND);
		repaint_region(ev, &repaint, &surface_blend);
		gs->used_in_output_repaint = true;
	}
#if HKMC_CCIC
#if QT_FAST_BLUR
	if (output->enable_fbo && ev->blur_radius > 0.0f && (ev->blur_end == 1 || ev->blur_gid < 1)) {
	#if DBG_BLUR_LOG
		weston_log("%s blur info : gid(%d), end(%d), radius(%.2f)\n",
			__func__, ev->blur_gid, ev->blur_end, ev->blur_radius);
	#endif
		repaint_blur(ev, output, &repaint);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
#elif QT_GAUSSIAN_BLUR
	if (output->enable_fbo && ev->blur_radius > 0.0f && (ev->blur_end == 1 || ev->blur_gid < 1)) {
		repaint_gaussian_blur(ev, output, &repaint, true);
		repaint_gaussian_blur(ev, output, &repaint, false);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
#endif
#endif
	pixman_region32_fini(&surface_blend);
	pixman_region32_fini(&surface_opaque);
out:
	pixman_region32_fini(&repaint);
	if (replaced_shader)
		gs->shader = replaced_shader;
}
#if HKMC_CCIC
#if defined(TILT_ENABLE) || defined(QT_FAST_BLUR) || defined(QT_GAUSSIAN_BLUR)
bool is_rear_screen(char *output_name) {
	if (strcmp(output_name, "HDMI-A-2") == 0) {
		return true;
	} else {
		return false;
	}
}
#endif
#if defined(QT_FAST_BLUR) || defined(QT_GAUSSIAN_BLUR)
void convert_tcoord(struct weston_output *output, pixman_box32_t *region, GLfloat *tcoord)
{
	int32_t width = output->current_mode->width;
	int32_t height = output->current_mode->height;
	if (region->x1 <= 0 || width <= 0) {
		tcoord[0] = 0;
	} else {
		tcoord[0] = (float)(region->x1) / width;
	}
	if (region->x2 <= 0 || width <= 0) {
		tcoord[1] = 0;
	} else {
		tcoord[1] = (float)(region->x2 + region->x1) / width;
	}
	int32_t y1_dividend = height - region->y1;
	if (y1_dividend <= 0 || height <= 0) {
		tcoord[3] = 0;
	} else {
		tcoord[3] = (float)y1_dividend / height;
	}
	int32_t y2_dividend = height - region->y2;
	if (y2_dividend <= 0 || height <= 0) {
		tcoord[2] = 0;
	} else {
		tcoord[2] = (float)y2_dividend / height;
	}
#if DBG_BLUR_LOG
	weston_log("region x %d,%d tcoord %.2f,%.2f,%.2f,%.2f\n",
		region->x1, region->x2, tcoord[0], tcoord[1], tcoord[2], tcoord[3]);
#endif
}
void scaling_blur_coord(struct weston_output *output, pixman_box32_t *rects, GLfloat *verts) {
	int32_t width = output->current_mode->width;
	int32_t height = output->current_mode->height;
	int32_t x1 = (*rects).x1 - output->x;
	int32_t x2 = (*rects).x2 - output->x;
	float vert_x1 = (float)x1 / (float)width;
	float vert_x2 = (float)x2 / (float)width;
	float vert_y1 = (float)(height - (*rects).y1) / (float)height;
	float vert_y2 = (float)(height - (*rects).y2) / (float)height;
	verts[0] = verts[6] = vert_x1;
	verts[1] = verts[3] = vert_y1;
	verts[2] = verts[4] = vert_x2;
	verts[5] = verts[7] = vert_y2;
}
void create_fbo(struct weston_output *output, struct gl_output_state *go)
{
	int i;
	int width = output->current_mode->width;
	int height = output->current_mode->height;
	for (i = 0; i < FBO_REPAINT_MAX; i++) {
		glGenFramebuffers(1, &go->blur_fbo[i]);
		glBindFramebuffer(GL_FRAMEBUFFER, go->blur_fbo[i]);
		glGenTextures(1, &go->blur_tex[i]);
		glBindTexture(GL_TEXTURE_2D, go->blur_tex[i]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_EXT, width, height, 0,
			GL_RGBA, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, go->blur_tex[i], 0);
		glBindTexture(GL_TEXTURE_2D, 0);
		int fb_status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		if (fb_status != GL_FRAMEBUFFER_COMPLETE) {
			weston_log("Unable to create blur FBO %d\n", i);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			return;
		}
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void
bind_fbo(struct weston_view *ev, struct weston_output *output)
{
	if (ev->blur_radius > 0.0f) {
	struct gl_output_state *go = get_output_state(output);
	glBindFramebuffer(GL_FRAMEBUFFER, go->blur_fbo[0]);
#if DBG_BLUR_LOG
	weston_log("%s : gid(%d), end(%d), radius(%.2f)\n",
		__func__, ev->blur_gid, ev->blur_end, ev->blur_radius);
#endif
	} else {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
}
#endif//end of if_defined
#if QT_FAST_BLUR
GLfloat get_weight(GLfloat v) {
	if (v <= 0.0) {
		return 1.0;
	} else if (v >= 0.5) {
		return 0.0;
	}
	return 1.0 - v * 2.0;
}
static void
repaint_fbo(struct weston_view *ev, struct weston_output *output,
	pixman_region32_t *damage, int step, float x_step, float y_step)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	struct gl_shader *shader = &gr->texture_shader_internal_blur;
	int32_t viewport_x = go->borders[GL_RENDERER_BORDER_LEFT].width;
	int32_t viewport_y = go->borders[GL_RENDERER_BORDER_BOTTOM].height;
	int32_t width = output->current_mode->width;
	int32_t height = output->current_mode->height;
	pixman_box32_t *rects;
	int n_rects;
	int i;
	static const GLfloat proj[16] = {
		2.0f,  0.0f, 0.0f, 0.0f,
		0.0f,  2.0f, 0.0f, 0.0f,
		0.0f,  0.0f, 1.0f, 0.0f,
		-1.0f, -1.0f, 0.0f, 1.0f
	};
	int current_fbo_id = go->blur_fbo[step];
	int before_texture_id = go->blur_tex[step - 1];
	use_shader(gr, shader);
	glUniform1f(shader->x_step_uniform, x_step);
	glUniform1f(shader->y_step_uniform, y_step);
	glUniform1i(shader->source_texture_uniform ,0);
	glUniformMatrix4fv(shader->proj_uniform, 1, GL_FALSE, proj);
	if (is_rear_screen(output->name)) {
		if (ev->blur_region.x1 >= REAR_EACH_WIDTH) { //Rear Right
			glUniform2i(shader->rear_screen_uniform, 1, 1);
#if DBG_BLUR_LOG
			weston_log("%s: Send Rear Right Screen Uniform\n", __func__);
#endif
		} else { //Rear Left
			glUniform2i(shader->rear_screen_uniform, 1, 0);
#if DBG_BLUR_LOG
			weston_log("%s: Send Rear Left Screen Uniform\n", __func__);
#endif
		}
	} else {
		glUniform2i(shader->rear_screen_uniform, 0, 0);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, current_fbo_id);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, before_texture_id);
	glViewport(viewport_x, viewport_y, width, height);
	rects = pixman_region32_rectangles(damage, &n_rects);
	GLfloat verts[4 * 2] = { 0.0f };
	for (i = 0; i < n_rects; i++) {
		scaling_blur_coord(output, rects, verts);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, verts);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, verts);
		glEnableVertexAttribArray(1);
		glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
#if DBG_BLUR_LOG
		weston_log("FBO render [%d %d] %.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
			(*rects).x1, (*rects).x2, verts[0], verts[1], verts[2], verts[3], verts[4], verts[5], verts[6], verts[7]);
#endif
	}
	glBindTexture(GL_TEXTURE_2D, 0);
}
static void
repaint_composition_fbo(struct weston_output *output, struct weston_view *ev, pixman_region32_t *damage)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	struct gl_shader *shader = &gr->texture_shader_fastblur;
	int32_t viewport_x = go->borders[GL_RENDERER_BORDER_LEFT].width;
	int32_t viewport_y = go->borders[GL_RENDERER_BORDER_BOTTOM].height;
	int32_t width = output->current_mode->width;
	int32_t height = output->current_mode->height;
	pixman_box32_t *rects;
	int n_rects;
	int i;
	static const GLfloat proj[16] = {
		2.0f,  0.0f, 0.0f, 0.0f,
		0.0f,  2.0f, 0.0f, 0.0f,
		0.0f,  0.0f, 1.0f, 0.0f,
		-1.0f, -1.0f, 0.0f, 1.0f
	};
	use_shader(gr, shader);
	shader_uniforms(shader, ev, output);
	GLfloat regions[4];
	convert_tcoord(output, &ev->blur_region, regions);
	glUniform4fv(shader->blur_region_uniform, 1, regions);
	glViewport(viewport_x, viewport_y, width, height);
	glUniform1i(shader->source_texture_uniform ,0);
	glUniformMatrix4fv(shader->proj_uniform, 1, GL_FALSE, proj);
	GLfloat weight[6] = {0.0, };
	GLfloat lod = sqrt(ev->blur_radius / 64.0) * 1.2 - 0.2;
	weight[0] = get_weight(fabs(lod - 0.100));
	weight[1] = get_weight(fabs(lod - 0.300));
	weight[2] = get_weight(fabs(lod - 0.500));
	weight[3] = get_weight(fabs(lod - 0.700));
	weight[4] = get_weight(fabs(lod - 0.900));
	GLfloat sum = weight[0] + weight[1] + weight[2] + weight[3] + weight[4];
	// 1 ~ 5 bind weight
	for (i = 0; i < 5; i++) {
		weight[i] /= sum;
		glUniform1f(shader->weight_uniform[i] , weight[i]);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	// 0 ~ 5 bind texture
	for (i = 0; i < FBO_REPAINT_MAX; i++) {
		glUniform1i(shader->source_uniform[i] ,i);
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(GL_TEXTURE_2D, go->blur_tex[i]);
	}
	glEnable(GL_BLEND);
	glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	rects = pixman_region32_rectangles(damage, &n_rects);
	GLfloat verts[4 * 2] = { 0.0f };
	scaling_blur_coord(output, rects, verts);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(1);
	glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
#if DBG_BLUR_LOG
	weston_log("Composition render %p [%d %d] %.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
		ev, (*rects).x1, (*rects).x2, verts[0], verts[1], verts[2], verts[3], verts[4], verts[5], verts[6], verts[7]);
#endif
	glDisable(GL_BLEND);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	glBindFramebuffer(GL_FRAMEBUFFER, go->blur_fbo[0]);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void repaint_blur(struct weston_view *ev, struct weston_output *output, pixman_region32_t *damage)
{
	// Render fbo Texture1~5, blur image.
	GLfloat blur_width = output->current_mode->width;
	GLfloat blur_height = output->current_mode->height;
	int step;
	for (step = 1; step < FBO_REPAINT_MAX; step++) {
		if (step == 1) {
			blur_width = ceil(blur_width / 32.0) * 32.0;
			blur_height = ceil(blur_height / 32.0) * 32.0;
		} else {
			blur_width /= 2;
			blur_height /= 2;
		}
		float x_step = 1.0f / blur_width;
		float y_step = 1.0f / blur_height;
		repaint_fbo(ev, output, damage, step, x_step, y_step);
	}
	// Composite Texutre 1~5.
	repaint_composition_fbo(output, ev, damage);
}
//end of QT_FAST_BLUR
#elif QT_GAUSSIAN_BLUR
void
set_gauss_properties(struct weston_view *ev, struct weston_output *output, struct gl_shader *shader,
			bool is_horizontal, float *spread, uint32_t *samples, float (* dirstep)[2])
{
	double deviation = ((int)ev->blur_radius + 1) / 3.3333;
	uint32_t requested_samples = (int)ev->blur_radius * 2 + 1;
	*samples = 1 + requested_samples / 2;
	uint32_t radius = requested_samples / 4;
	build_gauss_sample_points(shader, *samples, radius, deviation);
	uint32_t kernel_radius = *samples / 2;
	*spread = (float)radius / (float)kernel_radius;
	if (isnan(*spread) != 0) *spread = 0.0;
	int padded_tex_width = output->current_mode->width;
	int padded_tex_height = output->current_mode->height;
	float dpr = (float)padded_tex_width /(float)padded_tex_height;
	if(is_horizontal) {
		dirstep[0][0] = 1.0 / (padded_tex_width * dpr);
		dirstep[0][1] = 0;
	} else {
		dirstep[0][0] = 0;
		dirstep[0][1] = 1.0 / (padded_tex_height * dpr);
	}
}
void
repaint_gaussian_blur(struct weston_view *ev, struct weston_output *output,
			pixman_region32_t *damage, bool is_horizontal)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	struct gl_shader *shader = &gr->texture_shader_gaussianblur;
	int32_t viewport_x = go->borders[GL_RENDERER_BORDER_LEFT].width;
	int32_t viewport_y = go->borders[GL_RENDERER_BORDER_BOTTOM].height;
	int32_t width = output->current_mode->width;
	int32_t height = output->current_mode->height;
	pixman_box32_t *rects;
	int n_rects;
	int i;
	static const GLfloat proj[16] = {
		2.0f,  0.0f, 0.0f, 0.0f,
		0.0f,  2.0f, 0.0f, 0.0f,
		0.0f,  0.0f, 1.0f, 0.0f,
		-1.0f, -1.0f, 0.0f, 1.0f
	};
	use_shader(gr, shader);
	if(is_horizontal == false) {
		shader_uniforms(shader, ev, output);
	}
	GLfloat regions[4];
	convert_tcoord(output, &ev->blur_region, regions);
	glUniform4fv(shader->blur_region_uniform, 1, regions);
	float spread = 0.0;
	uint32_t samples = 0;
	float dirstep[2] = {0.0};
	set_gauss_properties(ev, output, shader, is_horizontal, &spread, &samples, &dirstep);
	glUniform1f(shader->spread_uniform, spread);
	glUniform2f(shader->dirstep_uniform, dirstep[0], dirstep[1]);
	glUniform1i(shader->samples_uniform, samples);
#if DBG_BLUR_LOG
	weston_log("%s(%d) spread(%.4f)  dirstep(%.4f/%.4f) samples(%d)\n",
		__func__, is_horizontal?0:1, spread, dirstep[0], dirstep[1], samples);
#endif
	glUniform1i(shader->source_texture_uniform, 0);
	glUniformMatrix4fv(shader->proj_uniform, 1, GL_FALSE, proj);
	//if (is_rear_screen(output->name)) //TODO: rear screen exception code
	if(is_horizontal == true) {
		glBindFramebuffer(GL_FRAMEBUFFER, go->blur_fbo[1]);
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, go->blur_tex[0]);
	} else {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, go->blur_tex[1]);
		glEnable(GL_BLEND);
		glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	}
	glViewport(viewport_x, viewport_y, width, height);
	rects = pixman_region32_rectangles(damage, &n_rects);
	GLfloat verts[4 * 2] = { 0.0f };
	scaling_blur_coord(output, rects, verts);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(1);
	glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
	glBindTexture(GL_TEXTURE_2D, 0);
	if(is_horizontal == false) {
		glDisable(GL_BLEND);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(0);
		for(i = 0; i < FBO_REPAINT_MAX; i++) {
			glBindFramebuffer(GL_FRAMEBUFFER, go->blur_fbo[i]);
			glClearColor(0.0, 0.0, 0.0, 0.0);
			glClear(GL_COLOR_BUFFER_BIT);
		}
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
}
#endif // end of QT_GAUSSIAN_BLUR
#endif // HKMC_CCIC
static void
repaint_views(struct weston_output *output, pixman_region32_t *damage)
{
	struct weston_compositor *compositor = output->compositor;
	struct weston_view *view;
	wl_list_for_each_reverse(view, &compositor->view_list, link)
		if (view->plane == &compositor->primary_plane)
			draw_view(view, output, damage);
}
static int
gl_renderer_create_fence_fd(struct weston_output *output);
/* Updates the release fences of surfaces that were used in the current output
 * repaint. Should only be used from gl_renderer_repaint_output, so that the
 * information in gl_surface_state.used_in_output_repaint is accurate.
 */
static void
update_buffer_release_fences(struct weston_compositor *compositor,
			     struct weston_output *output)
{
	struct weston_view *view;
	wl_list_for_each_reverse(view, &compositor->view_list, link) {
		struct gl_surface_state *gs;
		struct weston_buffer_release *buffer_release;
		int fence_fd;
		if (view->plane != &compositor->primary_plane)
			continue;
		gs = get_surface_state(view->surface);
		buffer_release = gs->buffer_release_ref.buffer_release;
		if (!gs->used_in_output_repaint || !buffer_release)
			continue;
		fence_fd = gl_renderer_create_fence_fd(output);
		/* If we have a buffer_release then it means we support fences,
		 * and we should be able to create the release fence. If we
		 * can't, something has gone horribly wrong, so disconnect the
		 * client.
		 */
		if (fence_fd == -1) {
			linux_explicit_synchronization_send_server_error(
				buffer_release->resource,
				"Failed to create release fence");
			fd_clear(&buffer_release->fence_fd);
			continue;
		}
		/* At the moment it is safe to just replace the fence_fd,
		 * discarding the previous one:
		 *
		 * 1. If the previous fence fd represents a sync fence from
		 *    a previous repaint cycle, that fence fd is now not
		 *    sufficient to provide the release guarantee and should
		 *    be replaced.
		 *
		 * 2. If the fence fd represents a sync fence from another
		 *    output in the same repaint cycle, it's fine to replace
		 *    it since we are rendering to all outputs using the same
		 *    EGL context, so a fence issued for a later output rendering
		 *    is guaranteed to signal after fences for previous output
		 *    renderings.
		 *
		 * Note that the above is only valid if the buffer_release
		 * fences only originate from the GL renderer, which guarantees
		 * a total order of operations and fences.  If we introduce
		 * fences from other sources (e.g., plane out-fences), we will
		 * need to merge fences instead.
		 */
		fd_update(&buffer_release->fence_fd, fence_fd);
	}
}
static void
draw_output_border_texture(struct gl_output_state *go,
			   enum gl_renderer_border_side side,
			   int32_t x, int32_t y,
			   int32_t width, int32_t height)
{
	struct gl_border_image *img = &go->borders[side];
	static GLushort indices [] = { 0, 1, 3, 3, 1, 2 };
	if (!img->data) {
		if (img->tex) {
			glDeleteTextures(1, &img->tex);
			img->tex = 0;
		}
		return;
	}
	if (!img->tex) {
		glGenTextures(1, &img->tex);
		glBindTexture(GL_TEXTURE_2D, img->tex);
		glTexParameteri(GL_TEXTURE_2D,
				GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D,
				GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D,
				GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D,
				GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	} else {
		glBindTexture(GL_TEXTURE_2D, img->tex);
	}
	if (go->border_status & (1 << side)) {
		glPixelStorei(GL_UNPACK_ROW_LENGTH_EXT, 0);
		glPixelStorei(GL_UNPACK_SKIP_PIXELS_EXT, 0);
		glPixelStorei(GL_UNPACK_SKIP_ROWS_EXT, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_BGRA_EXT,
			     img->tex_width, img->height, 0,
			     GL_BGRA_EXT, GL_UNSIGNED_BYTE, img->data);
	}
	GLfloat texcoord[] = {
		0.0f, 0.0f,
		(GLfloat)img->width / (GLfloat)img->tex_width, 0.0f,
		(GLfloat)img->width / (GLfloat)img->tex_width, 1.0f,
		0.0f, 1.0f,
	};
	GLfloat verts[] = {
		x, y,
		x + width, y,
		x + width, y + height,
		x, y + height
	};
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, texcoord);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, indices);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
}
static int
output_has_borders(struct weston_output *output)
{
	struct gl_output_state *go = get_output_state(output);
	return go->borders[GL_RENDERER_BORDER_TOP].data ||
	       go->borders[GL_RENDERER_BORDER_RIGHT].data ||
	       go->borders[GL_RENDERER_BORDER_BOTTOM].data ||
	       go->borders[GL_RENDERER_BORDER_LEFT].data;
}
static void
draw_output_borders(struct weston_output *output,
		    enum gl_border_status border_status)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	struct gl_shader *shader = &gr->texture_shader_rgba;
	struct gl_border_image *top, *bottom, *left, *right;
	struct weston_matrix matrix;
	int full_width, full_height;
	if (border_status == BORDER_STATUS_CLEAN)
		return; /* Clean. Nothing to do. */
	top = &go->borders[GL_RENDERER_BORDER_TOP];
	bottom = &go->borders[GL_RENDERER_BORDER_BOTTOM];
	left = &go->borders[GL_RENDERER_BORDER_LEFT];
	right = &go->borders[GL_RENDERER_BORDER_RIGHT];
	full_width = output->current_mode->width + left->width + right->width;
	full_height = output->current_mode->height + top->height + bottom->height;
	glDisable(GL_BLEND);
	use_shader(gr, shader);
	glViewport(0, 0, full_width, full_height);
	weston_matrix_init(&matrix);
	weston_matrix_translate(&matrix, -full_width/2.0, -full_height/2.0, 0);
	weston_matrix_scale(&matrix, 2.0/full_width, -2.0/full_height, 1);
	glUniformMatrix4fv(shader->proj_uniform, 1, GL_FALSE, matrix.d);
	glUniform1i(shader->tex_uniforms[0], 0);
	glUniform1f(shader->alpha_uniform, 1);
	glActiveTexture(GL_TEXTURE0);
	if (border_status & BORDER_TOP_DIRTY)
		draw_output_border_texture(go, GL_RENDERER_BORDER_TOP,
					   0, 0,
					   full_width, top->height);
	if (border_status & BORDER_LEFT_DIRTY)
		draw_output_border_texture(go, GL_RENDERER_BORDER_LEFT,
					   0, top->height,
					   left->width, output->current_mode->height);
	if (border_status & BORDER_RIGHT_DIRTY)
		draw_output_border_texture(go, GL_RENDERER_BORDER_RIGHT,
					   full_width - right->width, top->height,
					   right->width, output->current_mode->height);
	if (border_status & BORDER_BOTTOM_DIRTY)
		draw_output_border_texture(go, GL_RENDERER_BORDER_BOTTOM,
					   0, full_height - bottom->height,
					   full_width, bottom->height);
}
static void
output_get_border_damage(struct weston_output *output,
			 enum gl_border_status border_status,
			 pixman_region32_t *damage)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_border_image *top, *bottom, *left, *right;
	int full_width, full_height;
	if (border_status == BORDER_STATUS_CLEAN)
		return; /* Clean. Nothing to do. */
	top = &go->borders[GL_RENDERER_BORDER_TOP];
	bottom = &go->borders[GL_RENDERER_BORDER_BOTTOM];
	left = &go->borders[GL_RENDERER_BORDER_LEFT];
	right = &go->borders[GL_RENDERER_BORDER_RIGHT];
	full_width = output->current_mode->width + left->width + right->width;
	full_height = output->current_mode->height + top->height + bottom->height;
	if (border_status & BORDER_TOP_DIRTY)
		pixman_region32_union_rect(damage, damage,
					   0, 0,
					   full_width, top->height);
	if (border_status & BORDER_LEFT_DIRTY)
		pixman_region32_union_rect(damage, damage,
					   0, top->height,
					   left->width, output->current_mode->height);
	if (border_status & BORDER_RIGHT_DIRTY)
		pixman_region32_union_rect(damage, damage,
					   full_width - right->width, top->height,
					   right->width, output->current_mode->height);
	if (border_status & BORDER_BOTTOM_DIRTY)
		pixman_region32_union_rect(damage, damage,
					   0, full_height - bottom->height,
					   full_width, bottom->height);
}
static void
output_get_damage(struct weston_output *output,
		  pixman_region32_t *buffer_damage, uint32_t *border_damage)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	EGLint buffer_age = 0;
	EGLBoolean ret;
	int i;
	if (gr->has_egl_buffer_age) {
		ret = eglQuerySurface(gr->egl_display, go->egl_surface,
				      EGL_BUFFER_AGE_EXT, &buffer_age);
		if (ret == EGL_FALSE) {
			weston_log("buffer age query failed.\n");
			gl_renderer_print_egl_error_state();
		}
	}
	if (buffer_age == 0 || buffer_age - 1 > BUFFER_DAMAGE_COUNT) {
		pixman_region32_copy(buffer_damage, &output->region);
		*border_damage = BORDER_ALL_DIRTY;
	} else {
		for (i = 0; i < buffer_age - 1; i++)
			*border_damage |= go->border_damage[(go->buffer_damage_index + i) % BUFFER_DAMAGE_COUNT];
		if (*border_damage & BORDER_SIZE_CHANGED) {
			/* If we've had a resize, we have to do a full
			 * repaint. */
			*border_damage |= BORDER_ALL_DIRTY;
			pixman_region32_copy(buffer_damage, &output->region);
		} else {
			for (i = 0; i < buffer_age - 1; i++)
				pixman_region32_union(buffer_damage,
						      buffer_damage,
						      &go->buffer_damage[(go->buffer_damage_index + i) % BUFFER_DAMAGE_COUNT]);
		}
	}
}
static void
output_rotate_damage(struct weston_output *output,
		     pixman_region32_t *output_damage,
		     enum gl_border_status border_status)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	if (!gr->has_egl_buffer_age)
		return;
	go->buffer_damage_index += BUFFER_DAMAGE_COUNT - 1;
	go->buffer_damage_index %= BUFFER_DAMAGE_COUNT;
	pixman_region32_copy(&go->buffer_damage[go->buffer_damage_index], output_damage);
	go->border_damage[go->buffer_damage_index] = border_status;
}
/**
 * Given a region in Weston's (top-left-origin) global co-ordinate space,
 * translate it to the co-ordinate space used by GL for our output
 * rendering. This requires shifting it into output co-ordinate space:
 * translating for output offset within the global co-ordinate space,
 * multiplying by output scale to get buffer rather than logical size.
 *
 * Finally, if borders are drawn around the output, we translate the area
 * to account for the border region around the outside, and add any
 * damage if the borders have been redrawn.
 *
 * @param output The output whose co-ordinate space we are after
 * @param global_region The affected region in global co-ordinate space
 * @param[out] rects Y-inverted quads in {x,y,w,h} order; caller must free
 * @param[out] nrects Number of quads (4x number of co-ordinates)
 */
static void
pixman_region_to_egl_y_invert(struct weston_output *output,
			      struct pixman_region32 *global_region,
			      EGLint **rects,
			      EGLint *nrects)
{
	struct gl_output_state *go = get_output_state(output);
	pixman_region32_t transformed;
	struct pixman_box32 *box;
	int buffer_height;
	EGLint *d;
	int i;
	/* Translate from global to output co-ordinate space. */
	pixman_region32_init(&transformed);
	pixman_region32_copy(&transformed, global_region);
	pixman_region32_translate(&transformed, -output->x, -output->y);
	weston_transformed_region(output->width, output->height,
				  output->transform,
				  output->current_scale,
				  &transformed, &transformed);
	/* If we have borders drawn around the output, shift our output damage
	 * to account for borders being drawn around the outside, adding any
	 * damage resulting from borders being redrawn. */
	if (output_has_borders(output)) {
		pixman_region32_translate(&transformed,
					  go->borders[GL_RENDERER_BORDER_LEFT].width,
					  go->borders[GL_RENDERER_BORDER_TOP].height);
		output_get_border_damage(output, go->border_status,
					 &transformed);
	}
	/* Convert from a Pixman region into {x,y,w,h} quads, flipping in the
	 * Y axis to account for GL's lower-left-origin co-ordinate space. */
	box = pixman_region32_rectangles(&transformed, nrects);
	*rects = malloc(*nrects * 4 * sizeof(EGLint));
	buffer_height = go->borders[GL_RENDERER_BORDER_TOP].height +
			output->current_mode->height +
			go->borders[GL_RENDERER_BORDER_BOTTOM].height;
	d = *rects;
	for (i = 0; i < *nrects; ++i) {
		*d++ = box[i].x1;
		*d++ = buffer_height - box[i].y2;
		*d++ = box[i].x2 - box[i].x1;
		*d++ = box[i].y2 - box[i].y1;
	}
	pixman_region32_fini(&transformed);
}
/* NOTE: We now allow falling back to ARGB gl visuals when XRGB is
 * unavailable, so we're assuming the background has no transparency
 * and that everything with a blend, like drop shadows, will have something
 * opaque (like the background) drawn underneath it.
 *
 * Depending on the underlying hardware, violating that assumption could
 * result in seeing through to another display plane.
 */
static void
gl_renderer_repaint_output(struct weston_output *output,
			      pixman_region32_t *output_damage)
{
	struct gl_output_state *go = get_output_state(output);
	struct weston_compositor *compositor = output->compositor;
	struct gl_renderer *gr = get_renderer(compositor);
	EGLBoolean ret;
	static int errored;
	/* areas we've damaged since we last used this buffer */
	pixman_region32_t previous_damage;
	/* total area we need to repaint this time */
	pixman_region32_t total_damage;
	enum gl_border_status border_status = BORDER_STATUS_CLEAN;
	struct weston_view *view;
	if (use_output(output) < 0)
		return;
	/* Clear the used_in_output_repaint flag, so that we can properly track
	 * which surfaces were used in this output repaint. */
	wl_list_for_each_reverse(view, &compositor->view_list, link) {
		if (view->plane == &compositor->primary_plane) {
			struct gl_surface_state *gs =
				get_surface_state(view->surface);
			gs->used_in_output_repaint = false;
		}
	}
	if (go->begin_render_sync != EGL_NO_SYNC_KHR)
		gr->destroy_sync(gr->egl_display, go->begin_render_sync);
	if (go->end_render_sync != EGL_NO_SYNC_KHR)
		gr->destroy_sync(gr->egl_display, go->end_render_sync);
	go->begin_render_sync = create_render_sync(gr);
	/* Calculate the viewport */
	glViewport(go->borders[GL_RENDERER_BORDER_LEFT].width,
		   go->borders[GL_RENDERER_BORDER_BOTTOM].height,
		   output->current_mode->width,
		   output->current_mode->height);
	/* Calculate the global GL matrix */
	go->output_matrix = output->matrix;
	weston_matrix_translate(&go->output_matrix,
				-(output->current_mode->width / 2.0),
				-(output->current_mode->height / 2.0), 0);
	weston_matrix_scale(&go->output_matrix,
			    2.0 / output->current_mode->width,
			    -2.0 / output->current_mode->height, 1);
#if 1 // nvidia
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
#endif
	/* In fan debug mode, redraw everything to make sure that we clear any
	 * fans left over from previous draws on this buffer.
	 * This precludes the use of EGL_EXT_swap_buffers_with_damage and
	 * EGL_KHR_partial_update, since we damage the whole area. */
	if (gr->fan_debug) {
		pixman_region32_t undamaged;
		pixman_region32_init(&undamaged);
		pixman_region32_subtract(&undamaged, &output->region,
					 output_damage);
		gr->fan_debug = false;
		repaint_views(output, &undamaged);
		gr->fan_debug = true;
		pixman_region32_fini(&undamaged);
	}
	/* previous_damage covers regions damaged in previous paints since we
	 * last used this buffer */
	pixman_region32_init(&previous_damage);
	pixman_region32_init(&total_damage); /* total area to redraw */
	/* Update previous_damage using buffer_age (if available), and store
	 * current damaged region for future use. */
	output_get_damage(output, &previous_damage, &border_status);
	output_rotate_damage(output, output_damage, go->border_status);
	/* Redraw both areas which have changed since we last used this buffer,
	 * as well as the areas we now want to repaint, to make sure the
	 * buffer is up to date. */
	pixman_region32_union(&total_damage, &previous_damage, output_damage);
	border_status |= go->border_status;
	if (gr->has_egl_partial_update && !gr->fan_debug) {
		int n_egl_rects;
		EGLint *egl_rects;
		/* For partial_update, we need to pass the region which has
		 * changed since we last rendered into this specific buffer;
		 * this is total_damage. */
		pixman_region_to_egl_y_invert(output, &total_damage,
					      &egl_rects, &n_egl_rects);
		gr->set_damage_region(gr->egl_display, go->egl_surface,
				      egl_rects, n_egl_rects);
		free(egl_rects);
	}
	repaint_views(output, &total_damage);
	pixman_region32_fini(&total_damage);
	pixman_region32_fini(&previous_damage);
	draw_output_borders(output, border_status);
	wl_signal_emit(&output->frame_signal, output_damage);
	go->end_render_sync = create_render_sync(gr);
	if (gr->swap_buffers_with_damage && !gr->fan_debug) {
		int n_egl_rects;
		EGLint *egl_rects;
		/* For swap_buffers_with_damage, we need to pass the region
		 * which has changed since the previous SwapBuffers on this
		 * surface - this is output_damage. */
		pixman_region_to_egl_y_invert(output, output_damage,
					      &egl_rects, &n_egl_rects);
		ret = gr->swap_buffers_with_damage(gr->egl_display,
						   go->egl_surface,
						   egl_rects, n_egl_rects);
		free(egl_rects);
	} else {
		ret = eglSwapBuffers(gr->egl_display, go->egl_surface);
	}
	if (ret == EGL_FALSE && !errored) {
		errored = 1;
		weston_log("Failed in eglSwapBuffers.\n");
		gl_renderer_print_egl_error_state();
	}
	go->border_status = BORDER_STATUS_CLEAN;
	/* We have to submit the render sync objects after swap buffers, since
	 * the objects get assigned a valid sync file fd only after a gl flush.
	 */
	timeline_submit_render_sync(gr, compositor, output,
				    go->begin_render_sync,
				    TIMELINE_RENDER_POINT_TYPE_BEGIN);
	timeline_submit_render_sync(gr, compositor, output, go->end_render_sync,
				    TIMELINE_RENDER_POINT_TYPE_END);
	update_buffer_release_fences(compositor, output);
}
#if 1 // nvidia
static int
gl_renderer_output_stream_flip(struct weston_output *output,
			       void *flip_data)
{
#if defined(EGL_NV_stream_attrib) && defined(EGL_EXT_stream_acquire_mode)
	struct gl_output_state *go = get_output_state(output);
	struct weston_compositor *compositor = output->compositor;
	struct gl_renderer *gr = get_renderer(compositor);
	EGLAttrib acquire_attribs[3] = { EGL_NONE };
#ifdef EGL_NV_output_drm_flip_event
	if (gr->has_egl_output_drm_flip_event) {
		acquire_attribs[0] = EGL_DRM_FLIP_EVENT_DATA_NV;
		acquire_attribs[1] = (EGLAttrib)flip_data;
		acquire_attribs[2] = EGL_NONE;
	}
#endif
	if (go->egl_stream != EGL_NO_STREAM_KHR)
		if (gr->stream_consumer_acquire_attrib(
						gr->egl_display,
						go->egl_stream,
						acquire_attribs) != EGL_TRUE)
			return -1;
	return 0;
#else
	return -1;
#endif
}
static int
gl_renderer_output_stream_enable_auto_acquire(struct weston_output *output)
{
	struct gl_output_state *go = get_output_state(output);
	struct weston_compositor *compositor = output->compositor;
	struct gl_renderer *gr = get_renderer(compositor);
	if (go->egl_stream == EGL_NO_STREAM_KHR) {
		weston_log("Enable auto acquire failed; stream not created\n");
		return -1;
	}
	/* Re-enable auto-acquire mode */
	if (!(gr->stream_attrib(gr->egl_display, go->egl_stream,
				EGL_CONSUMER_AUTO_ACQUIRE_EXT, EGL_TRUE))) {
		weston_log("Enable auto acquire failed; unable to set stream attribute\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	return 0;
}
#endif // nvidia
static int
gl_renderer_read_pixels(struct weston_output *output,
			       pixman_format_code_t format, void *pixels,
			       uint32_t x, uint32_t y,
			       uint32_t width, uint32_t height)
{
	GLenum gl_format;
	struct gl_output_state *go = get_output_state(output);
	x += go->borders[GL_RENDERER_BORDER_LEFT].width;
	y += go->borders[GL_RENDERER_BORDER_BOTTOM].height;
	switch (format) {
	case PIXMAN_a8r8g8b8:
		gl_format = GL_BGRA_EXT;
		break;
	case PIXMAN_a8b8g8r8:
		gl_format = GL_RGBA;
		break;
	default:
		return -1;
	}
	if (use_output(output) < 0)
		return -1;
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(x, y, width, height, gl_format,
		     GL_UNSIGNED_BYTE, pixels);
	return 0;
}
static GLenum gl_format_from_internal(GLenum internal_format)
{
	switch (internal_format) {
	case GL_R8_EXT:
		return GL_RED_EXT;
	case GL_RG8_EXT:
		return GL_RG_EXT;
	default:
		return internal_format;
	}
}
static void
gl_renderer_flush_damage(struct weston_surface *surface)
{
	struct gl_renderer *gr = get_renderer(surface->compositor);
	struct gl_surface_state *gs = get_surface_state(surface);
	struct weston_buffer *buffer = gs->buffer_ref.buffer;
	struct weston_view *view;
	bool texture_used;
	pixman_box32_t *rectangles;
	uint8_t *data;
	int i, j, n;
	pixman_region32_union(&gs->texture_damage,
			      &gs->texture_damage, &surface->damage);
	if (!buffer)
		return;
	/* Avoid upload, if the texture won't be used this time.
	 * We still accumulate the damage in texture_damage, and
	 * hold the reference to the buffer, in case the surface
	 * migrates back to the primary plane.
	 */
	texture_used = false;
	wl_list_for_each(view, &surface->views, surface_link) {
		if (view->plane == &surface->compositor->primary_plane) {
			texture_used = true;
			break;
		}
	}
	if (!texture_used)
		return;
	if (!pixman_region32_not_empty(&gs->texture_damage) &&
	    !gs->needs_full_upload)
		goto done;
	data = wl_shm_buffer_get_data(buffer->shm_buffer);
	if (!gr->has_unpack_subimage) {
		wl_shm_buffer_begin_access(buffer->shm_buffer);
		for (j = 0; j < gs->num_textures; j++) {
			glBindTexture(GL_TEXTURE_2D, gs->textures[j]);
			glTexImage2D(GL_TEXTURE_2D, 0,
				     gs->gl_format[j],
				     gs->pitch / gs->hsub[j],
				     buffer->height / gs->vsub[j],
				     0,
				     gl_format_from_internal(gs->gl_format[j]),
				     gs->gl_pixel_type,
				     data + gs->offset[j]);
		}
		wl_shm_buffer_end_access(buffer->shm_buffer);
		goto done;
	}
	if (gs->needs_full_upload) {
		glPixelStorei(GL_UNPACK_SKIP_PIXELS_EXT, 0);
		glPixelStorei(GL_UNPACK_SKIP_ROWS_EXT, 0);
		wl_shm_buffer_begin_access(buffer->shm_buffer);
		for (j = 0; j < gs->num_textures; j++) {
			glBindTexture(GL_TEXTURE_2D, gs->textures[j]);
			glPixelStorei(GL_UNPACK_ROW_LENGTH_EXT,
				      gs->pitch / gs->hsub[j]);
			glTexImage2D(GL_TEXTURE_2D, 0,
				     gs->gl_format[j],
				     gs->pitch / gs->hsub[j],
				     buffer->height / gs->vsub[j],
				     0,
				     gl_format_from_internal(gs->gl_format[j]),
				     gs->gl_pixel_type,
				     data + gs->offset[j]);
		}
		wl_shm_buffer_end_access(buffer->shm_buffer);
		goto done;
	}
	rectangles = pixman_region32_rectangles(&gs->texture_damage, &n);
	wl_shm_buffer_begin_access(buffer->shm_buffer);
	for (i = 0; i < n; i++) {
		pixman_box32_t r;
		r = weston_surface_to_buffer_rect(surface, rectangles[i]);
		for (j = 0; j < gs->num_textures; j++) {
			glBindTexture(GL_TEXTURE_2D, gs->textures[j]);
			glPixelStorei(GL_UNPACK_ROW_LENGTH_EXT,
				      gs->pitch / gs->hsub[j]);
			glPixelStorei(GL_UNPACK_SKIP_PIXELS_EXT,
				      r.x1 / gs->hsub[j]);
			glPixelStorei(GL_UNPACK_SKIP_ROWS_EXT,
				      r.y1 / gs->hsub[j]);
			glTexSubImage2D(GL_TEXTURE_2D, 0,
					r.x1 / gs->hsub[j],
					r.y1 / gs->vsub[j],
					(r.x2 - r.x1) / gs->hsub[j],
					(r.y2 - r.y1) / gs->vsub[j],
					gl_format_from_internal(gs->gl_format[j]),
					gs->gl_pixel_type,
					data + gs->offset[j]);
		}
	}
	wl_shm_buffer_end_access(buffer->shm_buffer);
done:
	pixman_region32_fini(&gs->texture_damage);
	pixman_region32_init(&gs->texture_damage);
	gs->needs_full_upload = false;
	weston_buffer_reference(&gs->buffer_ref, NULL);
	weston_buffer_release_reference(&gs->buffer_release_ref, NULL);
}
static void
ensure_textures(struct gl_surface_state *gs, int num_textures)
{
	int i;
	if (num_textures <= gs->num_textures)
		return;
	for (i = gs->num_textures; i < num_textures; i++) {
		glGenTextures(1, &gs->textures[i]);
		glBindTexture(gs->target, gs->textures[i]);
		glTexParameteri(gs->target,
				GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(gs->target,
				GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	}
	gs->num_textures = num_textures;
	glBindTexture(gs->target, 0);
}
static void
gl_renderer_attach_shm(struct weston_surface *es, struct weston_buffer *buffer,
		       struct wl_shm_buffer *shm_buffer)
{
	struct weston_compositor *ec = es->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(es);
	GLenum gl_format[3] = {0, 0, 0};
	GLenum gl_pixel_type;
	int pitch;
	int num_planes;
	buffer->shm_buffer = shm_buffer;
	buffer->width = wl_shm_buffer_get_width(shm_buffer);
	buffer->height = wl_shm_buffer_get_height(shm_buffer);
	num_planes = 1;
	gs->offset[0] = 0;
	gs->hsub[0] = 1;
	gs->vsub[0] = 1;
	switch (wl_shm_buffer_get_format(shm_buffer)) {
	case WL_SHM_FORMAT_XRGB8888:
		gs->shader = &gr->texture_shader_rgbx;
		pitch = wl_shm_buffer_get_stride(shm_buffer) / 4;
		gl_format[0] = GL_BGRA_EXT;
		gl_pixel_type = GL_UNSIGNED_BYTE;
		es->is_opaque = true;
		break;
	case WL_SHM_FORMAT_ARGB8888:
		gs->shader = &gr->texture_shader_rgba;
		pitch = wl_shm_buffer_get_stride(shm_buffer) / 4;
		gl_format[0] = GL_BGRA_EXT;
		gl_pixel_type = GL_UNSIGNED_BYTE;
		es->is_opaque = false;
		break;
	case WL_SHM_FORMAT_RGB565:
		gs->shader = &gr->texture_shader_rgbx;
		pitch = wl_shm_buffer_get_stride(shm_buffer) / 2;
		gl_format[0] = GL_RGB;
		gl_pixel_type = GL_UNSIGNED_SHORT_5_6_5;
		es->is_opaque = true;
		break;
	case WL_SHM_FORMAT_YUV420:
		gs->shader = &gr->texture_shader_y_u_v;
		pitch = wl_shm_buffer_get_stride(shm_buffer);
		gl_pixel_type = GL_UNSIGNED_BYTE;
		num_planes = 3;
		gs->offset[1] = gs->offset[0] + (pitch / gs->hsub[0]) *
				(buffer->height / gs->vsub[0]);
		gs->hsub[1] = 2;
		gs->vsub[1] = 2;
		gs->offset[2] = gs->offset[1] + (pitch / gs->hsub[1]) *
				(buffer->height / gs->vsub[1]);
		gs->hsub[2] = 2;
		gs->vsub[2] = 2;
		if (gr->has_gl_texture_rg) {
			gl_format[0] = GL_R8_EXT;
			gl_format[1] = GL_R8_EXT;
			gl_format[2] = GL_R8_EXT;
		} else {
			gl_format[0] = GL_LUMINANCE;
			gl_format[1] = GL_LUMINANCE;
			gl_format[2] = GL_LUMINANCE;
		}
		es->is_opaque = true;
		break;
	case WL_SHM_FORMAT_NV12:
		pitch = wl_shm_buffer_get_stride(shm_buffer);
		gl_pixel_type = GL_UNSIGNED_BYTE;
		num_planes = 2;
		gs->offset[1] = gs->offset[0] + (pitch / gs->hsub[0]) *
				(buffer->height / gs->vsub[0]);
		gs->hsub[1] = 2;
		gs->vsub[1] = 2;
		if (gr->has_gl_texture_rg) {
			gs->shader = &gr->texture_shader_y_uv;
			gl_format[0] = GL_R8_EXT;
			gl_format[1] = GL_RG8_EXT;
		} else {
			gs->shader = &gr->texture_shader_y_xuxv;
			gl_format[0] = GL_LUMINANCE;
			gl_format[1] = GL_LUMINANCE_ALPHA;
		}
		es->is_opaque = true;
		break;
	case WL_SHM_FORMAT_YUYV:
		gs->shader = &gr->texture_shader_y_xuxv;
		pitch = wl_shm_buffer_get_stride(shm_buffer) / 2;
		gl_pixel_type = GL_UNSIGNED_BYTE;
		num_planes = 2;
		gs->offset[1] = 0;
		gs->hsub[1] = 2;
		gs->vsub[1] = 1;
		if (gr->has_gl_texture_rg)
			gl_format[0] = GL_RG8_EXT;
		else
			gl_format[0] = GL_LUMINANCE_ALPHA;
		gl_format[1] = GL_BGRA_EXT;
		es->is_opaque = true;
		break;
	default:
		weston_log("warning: unknown shm buffer format: %08x\n",
			   wl_shm_buffer_get_format(shm_buffer));
		return;
	}
	/* Only allocate a texture if it doesn't match existing one.
	 * If a switch from DRM allocated buffer to a SHM buffer is
	 * happening, we need to allocate a new texture buffer. */
	if (pitch != gs->pitch ||
	    buffer->height != gs->height ||
	    gl_format[0] != gs->gl_format[0] ||
	    gl_format[1] != gs->gl_format[1] ||
	    gl_format[2] != gs->gl_format[2] ||
	    gl_pixel_type != gs->gl_pixel_type ||
	    gs->buffer_type != BUFFER_TYPE_SHM) {
		gs->pitch = pitch;
		gs->height = buffer->height;
		gs->target = GL_TEXTURE_2D;
		gs->gl_format[0] = gl_format[0];
		gs->gl_format[1] = gl_format[1];
		gs->gl_format[2] = gl_format[2];
		gs->gl_pixel_type = gl_pixel_type;
		gs->buffer_type = BUFFER_TYPE_SHM;
		gs->needs_full_upload = true;
		gs->y_inverted = true;
		gs->direct_display = false;
		gs->surface = es;
		ensure_textures(gs, num_planes);
	}
}
static void
gl_renderer_attach_egl(struct weston_surface *es, struct weston_buffer *buffer,
		       uint32_t format)
{
	struct weston_compositor *ec = es->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(es);
	EGLint attribs[3];
	int i, num_planes;
	buffer->legacy_buffer = (struct wl_buffer *)buffer->resource;
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_WIDTH, &buffer->width);
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_HEIGHT, &buffer->height);
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_WAYLAND_Y_INVERTED_WL, &buffer->y_inverted);
	for (i = 0; i < gs->num_images; i++) {
		egl_image_unref(gs->images[i]);
		gs->images[i] = NULL;
	}
	gs->num_images = 0;
	gs->target = GL_TEXTURE_2D;
	es->is_opaque = false;
	switch (format) {
	case EGL_TEXTURE_RGB:
		es->is_opaque = true;
		/* fallthrough */
	case EGL_TEXTURE_RGBA:
	default:
		num_planes = 1;
		gs->shader = &gr->texture_shader_rgba;
		break;
	case EGL_TEXTURE_EXTERNAL_WL:
		num_planes = 1;
		gs->target = GL_TEXTURE_EXTERNAL_OES;
		gs->shader = &gr->texture_shader_egl_external;
		break;
	case EGL_TEXTURE_Y_UV_WL:
		num_planes = 2;
		gs->shader = &gr->texture_shader_y_uv;
		es->is_opaque = true;
		break;
	case EGL_TEXTURE_Y_U_V_WL:
		num_planes = 3;
		gs->shader = &gr->texture_shader_y_u_v;
		es->is_opaque = true;
		break;
	case EGL_TEXTURE_Y_XUXV_WL:
		num_planes = 2;
		gs->shader = &gr->texture_shader_y_xuxv;
		es->is_opaque = true;
		break;
	}
	ensure_textures(gs, num_planes);
	for (i = 0; i < num_planes; i++) {
		attribs[0] = EGL_WAYLAND_PLANE_WL;
		attribs[1] = i;
		attribs[2] = EGL_NONE;
		gs->images[i] = egl_image_create(gr,
						 EGL_WAYLAND_BUFFER_WL,
						 buffer->legacy_buffer,
						 attribs);
		if (!gs->images[i]) {
			weston_log("failed to create img for plane %d\n", i);
			continue;
		}
		gs->num_images++;
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(gs->target, gs->textures[i]);
		gr->image_target_texture_2d(gs->target,
					    gs->images[i]->image);
	}
	gs->pitch = buffer->width;
	gs->height = buffer->height;
	gs->buffer_type = BUFFER_TYPE_EGL;
	gs->y_inverted = buffer->y_inverted;
}
static void
gl_renderer_destroy_dmabuf(struct linux_dmabuf_buffer *dmabuf)
{
	struct dmabuf_image *image = linux_dmabuf_buffer_get_user_data(dmabuf);
	dmabuf_image_destroy(image);
}
static struct egl_image *
import_simple_dmabuf(struct gl_renderer *gr,
                     struct dmabuf_attributes *attributes)
{
	struct egl_image *image;
	EGLint attribs[50];
	int atti = 0;
	bool has_modifier;
	/* This requires the Mesa commit in
	 * Mesa 10.3 (08264e5dad4df448e7718e782ad9077902089a07) or
	 * Mesa 10.2.7 (55d28925e6109a4afd61f109e845a8a51bd17652).
	 * Otherwise Mesa closes the fd behind our back and re-importing
	 * will fail.
	 * https://bugs.freedesktop.org/show_bug.cgi?id=76188
	 */
	attribs[atti++] = EGL_WIDTH;
	attribs[atti++] = attributes->width;
	attribs[atti++] = EGL_HEIGHT;
	attribs[atti++] = attributes->height;
	attribs[atti++] = EGL_LINUX_DRM_FOURCC_EXT;
	attribs[atti++] = attributes->format;
	if (attributes->modifier[0] != DRM_FORMAT_MOD_INVALID) {
		if (!gr->has_dmabuf_import_modifiers)
			return NULL;
		has_modifier = true;
	} else {
		has_modifier = false;
	}
	if (attributes->n_planes > 0) {
		attribs[atti++] = EGL_DMA_BUF_PLANE0_FD_EXT;
		attribs[atti++] = attributes->fd[0];
		attribs[atti++] = EGL_DMA_BUF_PLANE0_OFFSET_EXT;
		attribs[atti++] = attributes->offset[0];
		attribs[atti++] = EGL_DMA_BUF_PLANE0_PITCH_EXT;
		attribs[atti++] = attributes->stride[0];
		if (has_modifier) {
			attribs[atti++] = EGL_DMA_BUF_PLANE0_MODIFIER_LO_EXT;
			attribs[atti++] = attributes->modifier[0] & 0xFFFFFFFF;
			attribs[atti++] = EGL_DMA_BUF_PLANE0_MODIFIER_HI_EXT;
			attribs[atti++] = attributes->modifier[0] >> 32;
		}
	}
	if (attributes->n_planes > 1) {
		attribs[atti++] = EGL_DMA_BUF_PLANE1_FD_EXT;
		attribs[atti++] = attributes->fd[1];
		attribs[atti++] = EGL_DMA_BUF_PLANE1_OFFSET_EXT;
		attribs[atti++] = attributes->offset[1];
		attribs[atti++] = EGL_DMA_BUF_PLANE1_PITCH_EXT;
		attribs[atti++] = attributes->stride[1];
		if (has_modifier) {
			attribs[atti++] = EGL_DMA_BUF_PLANE1_MODIFIER_LO_EXT;
			attribs[atti++] = attributes->modifier[1] & 0xFFFFFFFF;
			attribs[atti++] = EGL_DMA_BUF_PLANE1_MODIFIER_HI_EXT;
			attribs[atti++] = attributes->modifier[1] >> 32;
		}
	}
	if (attributes->n_planes > 2) {
		attribs[atti++] = EGL_DMA_BUF_PLANE2_FD_EXT;
		attribs[atti++] = attributes->fd[2];
		attribs[atti++] = EGL_DMA_BUF_PLANE2_OFFSET_EXT;
		attribs[atti++] = attributes->offset[2];
		attribs[atti++] = EGL_DMA_BUF_PLANE2_PITCH_EXT;
		attribs[atti++] = attributes->stride[2];
		if (has_modifier) {
			attribs[atti++] = EGL_DMA_BUF_PLANE2_MODIFIER_LO_EXT;
			attribs[atti++] = attributes->modifier[2] & 0xFFFFFFFF;
			attribs[atti++] = EGL_DMA_BUF_PLANE2_MODIFIER_HI_EXT;
			attribs[atti++] = attributes->modifier[2] >> 32;
		}
	}
	if (gr->has_dmabuf_import_modifiers) {
		if (attributes->n_planes > 3) {
			attribs[atti++] = EGL_DMA_BUF_PLANE3_FD_EXT;
			attribs[atti++] = attributes->fd[3];
			attribs[atti++] = EGL_DMA_BUF_PLANE3_OFFSET_EXT;
			attribs[atti++] = attributes->offset[3];
			attribs[atti++] = EGL_DMA_BUF_PLANE3_PITCH_EXT;
			attribs[atti++] = attributes->stride[3];
			attribs[atti++] = EGL_DMA_BUF_PLANE3_MODIFIER_LO_EXT;
			attribs[atti++] = attributes->modifier[3] & 0xFFFFFFFF;
			attribs[atti++] = EGL_DMA_BUF_PLANE3_MODIFIER_HI_EXT;
			attribs[atti++] = attributes->modifier[3] >> 32;
		}
	}
	attribs[atti++] = EGL_NONE;
	image = egl_image_create(gr, EGL_LINUX_DMA_BUF_EXT, NULL,
				 attribs);
	return image;
}
/* The kernel header drm_fourcc.h defines the DRM formats below.  We duplicate
 * some of the definitions here so that building Weston won't require
 * bleeding-edge kernel headers.
 */
#ifndef DRM_FORMAT_R8
#define DRM_FORMAT_R8            fourcc_code('R', '8', ' ', ' ') /* [7:0] R */
#endif
#ifndef DRM_FORMAT_GR88
#define DRM_FORMAT_GR88          fourcc_code('G', 'R', '8', '8') /* [15:0] G:R 8:8 little endian */
#endif
#ifndef DRM_FORMAT_XYUV8888
#define DRM_FORMAT_XYUV8888      fourcc_code('X', 'Y', 'U', 'V') /* [31:0] X:Y:Cb:Cr 8:8:8:8 little endian */
#endif
struct yuv_format_descriptor yuv_formats[] = {
	{
		.format = DRM_FORMAT_YUYV,
		.input_planes = 1,
		.output_planes = 2,
		.texture_type = TEXTURE_Y_XUXV_WL,
		{{
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_GR88,
			.plane_index = 0
		}, {
			.width_divisor = 2,
			.height_divisor = 1,
			.format = DRM_FORMAT_ARGB8888,
			.plane_index = 0
		}}
	}, {
		.format = DRM_FORMAT_NV12,
		.input_planes = 2,
		.output_planes = 2,
		.texture_type = TEXTURE_Y_UV_WL,
		{{
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_R8,
			.plane_index = 0
		}, {
			.width_divisor = 2,
			.height_divisor = 2,
			.format = DRM_FORMAT_GR88,
			.plane_index = 1
		}}
	}, {
		.format = DRM_FORMAT_YUV420,
		.input_planes = 3,
		.output_planes = 3,
		.texture_type = TEXTURE_Y_U_V_WL,
		{{
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_R8,
			.plane_index = 0
		}, {
			.width_divisor = 2,
			.height_divisor = 2,
			.format = DRM_FORMAT_R8,
			.plane_index = 1
		}, {
			.width_divisor = 2,
			.height_divisor = 2,
			.format = DRM_FORMAT_R8,
			.plane_index = 2
		}}
	}, {
		.format = DRM_FORMAT_YUV444,
		.input_planes = 3,
		.output_planes = 3,
		.texture_type = TEXTURE_Y_U_V_WL,
		{{
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_R8,
			.plane_index = 0
		}, {
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_R8,
			.plane_index = 1
		}, {
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_R8,
			.plane_index = 2
		}}
	}, {
		.format = DRM_FORMAT_XYUV8888,
		.input_planes = 1,
		.output_planes = 1,
		.texture_type = TEXTURE_XYUV_WL,
		{{
			.width_divisor = 1,
			.height_divisor = 1,
			.format = DRM_FORMAT_XBGR8888,
			.plane_index = 0
		}}
	}
};
static struct egl_image *
import_dmabuf_single_plane(struct gl_renderer *gr,
                           const struct dmabuf_attributes *attributes,
                           struct yuv_plane_descriptor *descriptor)
{
	struct dmabuf_attributes plane;
	struct egl_image *image;
	char fmt[4];
	plane.width = attributes->width / descriptor->width_divisor;
	plane.height = attributes->height / descriptor->height_divisor;
	plane.format = descriptor->format;
	plane.n_planes = 1;
	plane.fd[0] = attributes->fd[descriptor->plane_index];
	plane.offset[0] = attributes->offset[descriptor->plane_index];
	plane.stride[0] = attributes->stride[descriptor->plane_index];
	plane.modifier[0] = attributes->modifier[descriptor->plane_index];
	image = import_simple_dmabuf(gr, &plane);
	if (!image) {
		weston_log("Failed to import plane %d as %.4s\n",
		           descriptor->plane_index,
		           dump_format(descriptor->format, fmt));
		return NULL;
	}
	return image;
}
static bool
import_yuv_dmabuf(struct gl_renderer *gr,
                  struct dmabuf_image *image)
{
	unsigned i;
	int j;
	int ret;
	struct yuv_format_descriptor *format = NULL;
	struct dmabuf_attributes *attributes = &image->dmabuf->attributes;
	char fmt[4];
	for (i = 0; i < ARRAY_LENGTH(yuv_formats); ++i) {
		if (yuv_formats[i].format == attributes->format) {
			format = &yuv_formats[i];
			break;
		}
	}
	if (!format) {
		weston_log("Error during import, and no known conversion for format "
		           "%.4s in the renderer\n",
		           dump_format(attributes->format, fmt));
		return false;
	}
	if (attributes->n_planes != format->input_planes) {
		weston_log("%.4s dmabuf must contain %d plane%s (%d provided)\n",
		           dump_format(format->format, fmt),
		           format->input_planes,
		           (format->input_planes > 1) ? "s" : "",
		           attributes->n_planes);
		return false;
	}
	for (j = 0; j < format->output_planes; ++j) {
		image->images[j] = import_dmabuf_single_plane(gr, attributes,
		                                              &format->plane[j]);
		if (!image->images[j]) {
			while (j) {
				ret = egl_image_unref(image->images[--j]);
				assert(ret == 0);
			}
			return false;
		}
	}
	image->num_images = format->output_planes;
	switch (format->texture_type) {
	case TEXTURE_Y_XUXV_WL:
		image->shader = &gr->texture_shader_y_xuxv;
		break;
	case TEXTURE_Y_UV_WL:
		image->shader = &gr->texture_shader_y_uv;
		break;
	case TEXTURE_Y_U_V_WL:
		image->shader = &gr->texture_shader_y_u_v;
		break;
	case TEXTURE_XYUV_WL:
		image->shader = &gr->texture_shader_xyuv;
		break;
	default:
		assert(false);
	}
	return true;
}
static GLenum
choose_texture_target(struct dmabuf_attributes *attributes)
{
	if (attributes->n_planes > 1)
		return GL_TEXTURE_EXTERNAL_OES;
	switch (attributes->format & ~DRM_FORMAT_BIG_ENDIAN) {
	case DRM_FORMAT_YUYV:
	case DRM_FORMAT_YVYU:
	case DRM_FORMAT_UYVY:
	case DRM_FORMAT_VYUY:
	case DRM_FORMAT_AYUV:
	case DRM_FORMAT_XYUV8888:
		return GL_TEXTURE_EXTERNAL_OES;
	default:
		return GL_TEXTURE_2D;
	}
}
static struct dmabuf_image *
import_dmabuf(struct gl_renderer *gr,
	      struct linux_dmabuf_buffer *dmabuf)
{
	struct egl_image *egl_image;
	struct dmabuf_image *image;
	image = dmabuf_image_create();
	image->dmabuf = dmabuf;
	egl_image = import_simple_dmabuf(gr, &dmabuf->attributes);
	if (egl_image) {
		image->num_images = 1;
		image->images[0] = egl_image;
		image->import_type = IMPORT_TYPE_DIRECT;
		image->target = choose_texture_target(&dmabuf->attributes);
		switch (image->target) {
		case GL_TEXTURE_2D:
			image->shader = &gr->texture_shader_rgba;
			break;
		default:
			image->shader = &gr->texture_shader_egl_external;
		}
	} else {
		if (!import_yuv_dmabuf(gr, image)) {
			dmabuf_image_destroy(image);
			return NULL;
		}
		image->import_type = IMPORT_TYPE_GL_CONVERSION;
		image->target = GL_TEXTURE_2D;
	}
	return image;
}
static void
gl_renderer_query_dmabuf_formats(struct weston_compositor *wc,
				int **formats, int *num_formats)
{
	struct gl_renderer *gr = get_renderer(wc);
	static const int fallback_formats[] = {
		DRM_FORMAT_ARGB8888,
		DRM_FORMAT_XRGB8888,
		DRM_FORMAT_YUYV,
		DRM_FORMAT_NV12,
		DRM_FORMAT_YUV420,
		DRM_FORMAT_YUV444,
		DRM_FORMAT_XYUV8888,
	};
	bool fallback = false;
	EGLint num;
	assert(gr->has_dmabuf_import);
	if (!gr->has_dmabuf_import_modifiers ||
	    !gr->query_dmabuf_formats(gr->egl_display, 0, NULL, &num)) {
		num = gr->has_gl_texture_rg ? ARRAY_LENGTH(fallback_formats) : 2;
		fallback = true;
	}
	*formats = calloc(num, sizeof(int));
	if (*formats == NULL) {
		*num_formats = 0;
		return;
	}
	if (fallback) {
		memcpy(*formats, fallback_formats, num * sizeof(int));
		*num_formats = num;
		return;
	}
	if (!gr->query_dmabuf_formats(gr->egl_display, num, *formats, &num)) {
		*num_formats = 0;
		free(*formats);
		return;
	}
	*num_formats = num;
}
static void
gl_renderer_query_dmabuf_modifiers(struct weston_compositor *wc, int format,
					uint64_t **modifiers,
					int *num_modifiers)
{
	struct gl_renderer *gr = get_renderer(wc);
	int num;
	assert(gr->has_dmabuf_import);
	if (!gr->has_dmabuf_import_modifiers ||
		!gr->query_dmabuf_modifiers(gr->egl_display, format, 0, NULL,
					    NULL, &num) ||
		num == 0) {
		*num_modifiers = 0;
		return;
	}
	*modifiers = calloc(num, sizeof(uint64_t));
	if (*modifiers == NULL) {
		*num_modifiers = 0;
		return;
	}
	if (!gr->query_dmabuf_modifiers(gr->egl_display, format,
				num, *modifiers, NULL, &num)) {
		*num_modifiers = 0;
		free(*modifiers);
		return;
	}
	*num_modifiers = num;
}
static bool
gl_renderer_import_dmabuf(struct weston_compositor *ec,
			  struct linux_dmabuf_buffer *dmabuf)
{
	struct gl_renderer *gr = get_renderer(ec);
	struct dmabuf_image *image;
	int i;
	assert(gr->has_dmabuf_import);
	for (i = 0; i < dmabuf->attributes.n_planes; i++) {
		/* return if EGL doesn't support import modifiers */
		if (dmabuf->attributes.modifier[i] != DRM_FORMAT_MOD_INVALID)
			if (!gr->has_dmabuf_import_modifiers)
				return false;
		/* return if modifiers passed are unequal */
		if (dmabuf->attributes.modifier[i] !=
		    dmabuf->attributes.modifier[0])
			return false;
	}
	/* reject all flags we do not recognize or handle */
	if (dmabuf->attributes.flags & ~ZWP_LINUX_BUFFER_PARAMS_V1_FLAGS_Y_INVERT)
		return false;
	image = import_dmabuf(gr, dmabuf);
	if (!image)
		return false;
	wl_list_insert(&gr->dmabuf_images, &image->link);
	linux_dmabuf_buffer_set_user_data(dmabuf, image,
		gl_renderer_destroy_dmabuf);
	return true;
}
static bool
import_known_dmabuf(struct gl_renderer *gr,
                    struct dmabuf_image *image)
{
	switch (image->import_type) {
	case IMPORT_TYPE_DIRECT:
		image->images[0] = import_simple_dmabuf(gr, &image->dmabuf->attributes);
		if (!image->images[0])
			return false;
		image->num_images = 1;
		break;
	case IMPORT_TYPE_GL_CONVERSION:
		if (!import_yuv_dmabuf(gr, image))
			return false;
		break;
	default:
		weston_log("Invalid import type for dmabuf\n");
		return false;
	}
	return true;
}
static bool
dmabuf_is_opaque(struct linux_dmabuf_buffer *dmabuf)
{
	const struct pixel_format_info *info;
	info = pixel_format_get_info(dmabuf->attributes.format &
				     ~DRM_FORMAT_BIG_ENDIAN);
	if (!info)
		return false;
	return pixel_format_is_opaque(info);
}
static void
gl_renderer_attach_dmabuf(struct weston_surface *surface,
			  struct weston_buffer *buffer,
			  struct linux_dmabuf_buffer *dmabuf)
{
	struct gl_renderer *gr = get_renderer(surface->compositor);
	struct gl_surface_state *gs = get_surface_state(surface);
	struct dmabuf_image *image;
	int i;
	int ret;
	if (!gr->has_dmabuf_import) {
		linux_dmabuf_buffer_send_server_error(dmabuf,
				"EGL dmabuf import not supported");
		return;
	}
	buffer->width = dmabuf->attributes.width;
	buffer->height = dmabuf->attributes.height;
	/*
	 * GL-renderer uses the OpenGL convention of texture coordinates, where
	 * the origin is at bottom-left. Because dmabuf buffers have the origin
	 * at top-left, we must invert the Y_INVERT flag to get the image right.
	 */
	buffer->y_inverted =
		!(dmabuf->attributes.flags & ZWP_LINUX_BUFFER_PARAMS_V1_FLAGS_Y_INVERT);
	for (i = 0; i < gs->num_images; i++)
		egl_image_unref(gs->images[i]);
	gs->num_images = 0;
	gs->pitch = buffer->width;
	gs->height = buffer->height;
	gs->buffer_type = BUFFER_TYPE_EGL;
	gs->y_inverted = buffer->y_inverted;
	gs->direct_display = dmabuf->direct_display;
	surface->is_opaque = dmabuf_is_opaque(dmabuf);
	/*
	 * We try to always hold an imported EGLImage from the dmabuf
	 * to prevent the client from preventing re-imports. But, we also
	 * need to re-import every time the contents may change because
	 * GL driver's caching may need flushing.
	 *
	 * Here we release the cache reference which has to be final.
	 */
	if (dmabuf->direct_display)
		return;
	image = linux_dmabuf_buffer_get_user_data(dmabuf);
	/* The dmabuf_image should have been created during the import */
	assert(image != NULL);
	for (i = 0; i < image->num_images; ++i) {
		ret = egl_image_unref(image->images[i]);
		assert(ret == 0);
	}
	if (!import_known_dmabuf(gr, image)) {
		linux_dmabuf_buffer_send_server_error(dmabuf, "EGL dmabuf import failed");
		return;
	}
	gs->num_images = image->num_images;
	for (i = 0; i < gs->num_images; ++i)
		gs->images[i] = egl_image_ref(image->images[i]);
	gs->target = image->target;
	ensure_textures(gs, gs->num_images);
	for (i = 0; i < gs->num_images; ++i) {
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(gs->target, gs->textures[i]);
		gr->image_target_texture_2d(gs->target, gs->images[i]->image);
	}
	gs->shader = image->shader;
}
#if 1 // nvidia
static bool
gl_renderer_set_output_gamma(struct weston_output *output,
			     float r, float g, float b)
{
	EGLint err;
	struct gl_output_state *go = get_output_state(output);
	if (r < 0.0 || r > 5.0 ||
	    g < 0.0 || g > 5.0 ||
	    b < 0.0 || b > 5.0)
		return EGL_BAD_PARAMETER;
	go->gamma[0] = r;
	go->gamma[1] = g;
	go->gamma[2] = b;
	return EGL_SUCCESS;
}
static bool
gl_renderer_set_output_shade(struct weston_output *output,
			     float r, float g, float b)
{
	EGLint err;
	struct gl_output_state *go = get_output_state(output);
	if (r < 0.1 || r > 2.0 ||
	    g < 0.1 || g > 2.0 ||
	    b < 0.1 || b > 2.0)
		return EGL_BAD_PARAMETER;
	go->shade[0] = r;
	go->shade[1] = g;
	go->shade[2] = b;
	return EGL_SUCCESS;
}
static EGLint
gl_renderer_attach_eglstream_consumer_common(struct weston_surface *es,
					     struct wl_resource *wl_eglstream)
{
#ifdef EGL_NV_stream_attrib
	struct weston_compositor *ec = es->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(es);
	EGLStreamKHR stream = EGL_NO_STREAM_KHR;
	EGLAttrib stream_attribs[] = {
#ifdef EGL_WL_wayland_eglstream
		EGL_WAYLAND_EGLSTREAM_WL, (EGLAttrib)wl_eglstream,
#endif
		EGL_NONE
	};
	/* Check for required extensions. If they arent supported, there's no
	 *  way the given resource corresponds to an EGLStream */
	if (!gr->has_egl_stream_attrib ||
	    !gr->has_egl_stream_consumer_gltexture ||
	    !gr->has_egl_wayland_eglstream)
		return EGL_BAD_ACCESS;
	stream = gr->create_stream_attrib(gr->egl_display, stream_attribs);
	if (stream == EGL_NO_STREAM_KHR)
		return eglGetError();
	if (gs->egl_stream != EGL_NO_STREAM_KHR)
		gr->destroy_stream(gr->egl_display, gs->egl_stream);
	gs->egl_stream = stream;
	gs->shader = &gr->texture_shader_egl_external;
	gs->target = GL_TEXTURE_EXTERNAL_OES;
	glActiveTexture(GL_TEXTURE0);
	ensure_textures(gs, 2);
	glBindTexture(gs->target, gs->textures[1]);
	gs->new_stream = (gr->stream_consumer_gltexture(
					gr->egl_display,
					gs->egl_stream) == EGL_TRUE);
	if (!gs->new_stream) {
		EGLint err = eglGetError();
		weston_log("failed to set stream consumer\n");
		weston_log("EGL error: 0x%04lx\n", (long)err);
		gr->destroy_stream(gr->egl_display, gs->egl_stream);
		gs->egl_stream = EGL_NO_STREAM_KHR;
		return err;
	}
	weston_surface_set_realtime(gs->surface, true);
	return EGL_SUCCESS;
#else
	return EGL_BAD_ACCESS;
#endif
}
static bool
gl_renderer_attach_eglstream_consumer(struct weston_surface *es,
				      struct wl_resource *stream)
{
	EGLint err = gl_renderer_attach_eglstream_consumer_common(es, stream);
	return (err == EGL_SUCCESS);
}
/*
 * gl_renderer_attach_stream_texture
 *
 * Try to bind given <buffer> to an EGLStream. If the given buffer was already
 * bound, it will acquire next frame on the stream.
 *
 * Return true if the given <buffer> corresponds to an EGLStream; otherwise,
 * return false (if might be another kind of buffer).
 */
static bool
gl_renderer_attach_stream_texture(struct weston_surface *es,
				  struct weston_buffer *buffer)
{
#ifdef EGL_NV_stream_attrib
	struct weston_compositor *ec = es->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(es);
	EGLint stream_state = EGL_STREAM_STATE_EMPTY_KHR;
	EGLint err;
	err = gl_renderer_attach_eglstream_consumer_common(es, buffer->resource);
	switch (err) {
	case EGL_BAD_ACCESS:
		/* EGL_BAD_ACCESS is generated whenever buffer->resource
		 * does not corresponds to a stream */
		return false;
	case EGL_BAD_STREAM_KHR:
		/* EGL_BAD_STREAM_KHR is generated whenever
		 * buffer->resource corresponds to a previously created
		 * stream so we must have a valid stream handle already
		 * we can use to acquire next frame */
		break;
	case EGL_SUCCESS:
		/* EGL_SUCCESS is returned when the server stream endpoint is
		 * created and a consumer is attached successfully. This path is
		 * hit with old drivers which use wl_surface::attach + commit to
		 * request stream creation. */
		return true; /* buffer->resource is EGLStream */
	default:
		/* An unknown error was generated */
		assert(0);
		return false;
	}
	/* At this point we should have a valid stream handle */
	assert(gs->egl_stream != EGL_NO_STREAM_KHR);
	/* Check whether there are new frames available */
	if (gr->query_stream(gr->egl_display,
			     gs->egl_stream,
			     EGL_STREAM_STATE_KHR,
			     &stream_state) != EGL_TRUE) {
		weston_log("failed to query stream state\n");
		gl_renderer_print_egl_error_state();
		return true; /* buffer->resource is EGLStream */
	}
	/* If no new frame available, re-use last one */
	if (stream_state != EGL_STREAM_STATE_NEW_FRAME_AVAILABLE_KHR) {
		/* Fake size of last frame */
		buffer->width = gs->pitch;
		buffer->height = gs->height;
		return true; /* buffer->resource is EGLStream */
	}
	if (gr->stream_consumer_acquire_attrib(gr->egl_display,
					       gs->egl_stream,
					       NULL) != EGL_TRUE) {
		weston_log("failed to acquire buffer\n");
		gl_renderer_print_egl_error_state();
		return true; /* buffer->resource is EGLStream */
	}
	/* Swap textures if new stream was created */
	if (gs->new_stream) {
		GLuint tmp = gs->textures[0];
		gs->textures[0] = gs->textures[1];
		gs->textures[1] = tmp;
		gs->new_stream = false;
	}
	/* Update buffer and surface data */
	buffer->legacy_buffer = (void *)buffer->resource;
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_WIDTH, &buffer->width);
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_HEIGHT, &buffer->height);
	gr->query_buffer(gr->egl_display, buffer->legacy_buffer,
			 EGL_WAYLAND_Y_INVERTED_WL, &buffer->y_inverted);
	gs->pitch = buffer->width;
	gs->height = buffer->height;
	gs->buffer_type = BUFFER_TYPE_EGL;
	gs->y_inverted = buffer->y_inverted;
	return true; /* buffer->resource is EGLStream */
#else
	return false;
#endif
}
#endif // nvidia
static void
gl_renderer_attach(struct weston_surface *es, struct weston_buffer *buffer)
{
	struct weston_compositor *ec = es->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	struct gl_surface_state *gs = get_surface_state(es);
	struct wl_shm_buffer *shm_buffer;
	struct linux_dmabuf_buffer *dmabuf;
	EGLint format;
	int i;
#if 1 // nvidia
	bool is_realtime = false;
#endif
	weston_buffer_reference(&gs->buffer_ref, buffer);
	weston_buffer_release_reference(&gs->buffer_release_ref,
					es->buffer_release_ref.buffer_release);
	if (!buffer) {
		for (i = 0; i < gs->num_images; i++) {
			egl_image_unref(gs->images[i]);
			gs->images[i] = NULL;
		}
		gs->num_images = 0;
		glDeleteTextures(gs->num_textures, gs->textures);
		gs->num_textures = 0;
		gs->buffer_type = BUFFER_TYPE_NULL;
		gs->y_inverted = true;
		gs->direct_display = false;
		es->is_opaque = false;
#if 1 // nvidia
		if (gs->egl_stream != EGL_NO_STREAM_KHR) {
			gr->destroy_stream(gr->egl_display, gs->egl_stream);
			gs->egl_stream = EGL_NO_STREAM_KHR;
		}
		weston_surface_set_realtime(es, false);
#endif
		return;
	}
	shm_buffer = wl_shm_buffer_get(buffer->resource);
	if (shm_buffer)
		gl_renderer_attach_shm(es, buffer, shm_buffer);
#if 1 // nvidia
	else if (gl_renderer_attach_stream_texture(es, buffer)) {
		/* The stream texture is attached. The surface is marked as
		 * realtime so its frame callbacks are always called. */
		 is_realtime = true;
	}
#endif
	else if (gr->has_bind_display &&
		 gr->query_buffer(gr->egl_display, (void *)buffer->resource,
				  EGL_TEXTURE_FORMAT, &format))
		gl_renderer_attach_egl(es, buffer, format);
	else if ((dmabuf = linux_dmabuf_buffer_get(buffer->resource)))
		gl_renderer_attach_dmabuf(es, buffer, dmabuf);
	else {
		weston_log("unhandled buffer type!\n");
		if (gr->has_bind_display) {
		  weston_log("eglQueryWaylandBufferWL failed\n");
		  gl_renderer_print_egl_error_state();
		}
		weston_buffer_reference(&gs->buffer_ref, NULL);
		weston_buffer_release_reference(&gs->buffer_release_ref, NULL);
		gs->buffer_type = BUFFER_TYPE_NULL;
		gs->y_inverted = true;
		es->is_opaque = false;
                weston_buffer_send_server_error(buffer,
			"disconnecting due to unhandled buffer type");
	}
#if 1 // nvidia
	weston_surface_set_realtime(es, is_realtime);
#endif
}
static void
gl_renderer_surface_set_color(struct weston_surface *surface,
		 float red, float green, float blue, float alpha)
{
	struct gl_surface_state *gs = get_surface_state(surface);
	struct gl_renderer *gr = get_renderer(surface->compositor);
	gs->color[0] = red;
	gs->color[1] = green;
	gs->color[2] = blue;
	gs->color[3] = alpha;
	gs->buffer_type = BUFFER_TYPE_SOLID;
	gs->pitch = 1;
	gs->height = 1;
	gs->shader = &gr->solid_shader;
}
static void
gl_renderer_surface_get_content_size(struct weston_surface *surface,
				     int *width, int *height)
{
	struct gl_surface_state *gs = get_surface_state(surface);
	if (gs->buffer_type == BUFFER_TYPE_NULL) {
		*width = 0;
		*height = 0;
	} else {
		*width = gs->pitch;
		*height = gs->height;
	}
}
static uint32_t
pack_color(pixman_format_code_t format, float *c)
{
	uint8_t r = round(c[0] * 255.0f);
	uint8_t g = round(c[1] * 255.0f);
	uint8_t b = round(c[2] * 255.0f);
	uint8_t a = round(c[3] * 255.0f);
	switch (format) {
	case PIXMAN_a8b8g8r8:
		return (a << 24) | (b << 16) | (g << 8) | r;
	default:
		assert(0);
		return 0;
	}
}
static int
gl_renderer_surface_copy_content(struct weston_surface *surface,
				 void *target, size_t size,
				 int src_x, int src_y,
				 int width, int height)
{
	static const GLfloat verts[4 * 2] = {
		0.0f, 0.0f,
		1.0f, 0.0f,
		1.0f, 1.0f,
		0.0f, 1.0f
	};
	static const GLfloat projmat_normal[16] = { /* transpose */
		 2.0f,  0.0f, 0.0f, 0.0f,
		 0.0f,  2.0f, 0.0f, 0.0f,
		 0.0f,  0.0f, 1.0f, 0.0f,
		-1.0f, -1.0f, 0.0f, 1.0f
	};
	static const GLfloat projmat_yinvert[16] = { /* transpose */
		 2.0f,  0.0f, 0.0f, 0.0f,
		 0.0f, -2.0f, 0.0f, 0.0f,
		 0.0f,  0.0f, 1.0f, 0.0f,
		-1.0f,  1.0f, 0.0f, 1.0f
	};
	const pixman_format_code_t format = PIXMAN_a8b8g8r8;
	const size_t bytespp = 4; /* PIXMAN_a8b8g8r8 */
	const GLenum gl_format = GL_RGBA; /* PIXMAN_a8b8g8r8 little-endian */
	struct gl_renderer *gr = get_renderer(surface->compositor);
	struct gl_surface_state *gs = get_surface_state(surface);
	int cw, ch;
	GLuint fbo;
	GLuint tex;
	GLenum status;
	const GLfloat *proj;
	int i;
	gl_renderer_surface_get_content_size(surface, &cw, &ch);
	switch (gs->buffer_type) {
	case BUFFER_TYPE_NULL:
		return -1;
	case BUFFER_TYPE_SOLID:
		*(uint32_t *)target = pack_color(format, gs->color);
		return 0;
	case BUFFER_TYPE_SHM:
		gl_renderer_flush_damage(surface);
		/* fall through */
	case BUFFER_TYPE_EGL:
		break;
	}
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, cw, ch,
		     0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindTexture(GL_TEXTURE_2D, 0);
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
			       GL_TEXTURE_2D, tex, 0);
	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE) {
		weston_log("%s: fbo error: %#x\n", __func__, status);
		glDeleteFramebuffers(1, &fbo);
		glDeleteTextures(1, &tex);
		return -1;
	}
	glViewport(0, 0, cw, ch);
	glDisable(GL_BLEND);
	use_shader(gr, gs->shader);
	if (gs->y_inverted)
		proj = projmat_normal;
	else
		proj = projmat_yinvert;
	glUniformMatrix4fv(gs->shader->proj_uniform, 1, GL_FALSE, proj);
	glUniform1f(gs->shader->alpha_uniform, 1.0f);
	for (i = 0; i < gs->num_textures; i++) {
		glUniform1i(gs->shader->tex_uniforms[i], i);
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(gs->target, gs->textures[i]);
		glTexParameteri(gs->target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(gs->target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	}
	/* position: */
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(0);
	/* texcoord: */
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, verts);
	glEnableVertexAttribArray(1);
	glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	glPixelStorei(GL_PACK_ALIGNMENT, bytespp);
	glReadPixels(src_x, src_y, width, height, gl_format,
		     GL_UNSIGNED_BYTE, target);
	glDeleteFramebuffers(1, &fbo);
	glDeleteTextures(1, &tex);
	return 0;
}
static void
surface_state_destroy(struct gl_surface_state *gs, struct gl_renderer *gr)
{
	int i;
	wl_list_remove(&gs->surface_destroy_listener.link);
	wl_list_remove(&gs->renderer_destroy_listener.link);
	gs->surface->renderer_state = NULL;
	glDeleteTextures(gs->num_textures, gs->textures);
	for (i = 0; i < gs->num_images; i++)
		egl_image_unref(gs->images[i]);
	weston_buffer_reference(&gs->buffer_ref, NULL);
	weston_buffer_release_reference(&gs->buffer_release_ref, NULL);
	pixman_region32_fini(&gs->texture_damage);
#if 1 // nvidia
	if (gs->egl_stream != EGL_NO_STREAM_KHR)
		gr->destroy_stream(gr->egl_display, gs->egl_stream);
#endif
	free(gs);
}
static void
surface_state_handle_surface_destroy(struct wl_listener *listener, void *data)
{
	struct gl_surface_state *gs;
	struct gl_renderer *gr;
	gs = container_of(listener, struct gl_surface_state,
			  surface_destroy_listener);
	gr = get_renderer(gs->surface->compositor);
	surface_state_destroy(gs, gr);
}
static void
surface_state_handle_renderer_destroy(struct wl_listener *listener, void *data)
{
	struct gl_surface_state *gs;
	struct gl_renderer *gr;
	gr = data;
	gs = container_of(listener, struct gl_surface_state,
			  renderer_destroy_listener);
	surface_state_destroy(gs, gr);
}
static int
gl_renderer_create_surface(struct weston_surface *surface)
{
	struct gl_surface_state *gs;
	struct gl_renderer *gr = get_renderer(surface->compositor);
	gs = zalloc(sizeof *gs);
	if (gs == NULL)
		return -1;
	/* A buffer is never attached to solid color surfaces, yet
	 * they still go through texcoord computations. Do not divide
	 * by zero there.
	 */
	gs->pitch = 1;
	gs->y_inverted = true;
	gs->direct_display = false;
	gs->surface = surface;
#if 1 // nvidia
	gs->egl_stream = EGL_NO_STREAM_KHR;
#endif
	pixman_region32_init(&gs->texture_damage);
	surface->renderer_state = gs;
	gs->surface_destroy_listener.notify =
		surface_state_handle_surface_destroy;
	wl_signal_add(&surface->destroy_signal,
		      &gs->surface_destroy_listener);
	gs->renderer_destroy_listener.notify =
		surface_state_handle_renderer_destroy;
	wl_signal_add(&gr->destroy_signal,
		      &gs->renderer_destroy_listener);
	if (surface->buffer_ref.buffer) {
		gl_renderer_attach(surface, surface->buffer_ref.buffer);
		gl_renderer_flush_damage(surface);
	}
	return 0;
}
#if HKMC_CCIC
static const char vertex_shader[] =
	"uniform mat4 proj;\n"
	"uniform mat4 tilt;\n"
	"attribute vec2 position;\n"
	"attribute vec2 texcoord;\n"
	"varying vec2 v_texcoord;\n"
	"void main()\n"
	"{\n"
	"   gl_Position = tilt * proj * vec4(position, 0.0, 1.0);\n"
	"   v_texcoord = texcoord;\n"
	"}\n";
#if QT_FAST_BLUR
static const char texture_vertex_shader_internal_blur[] =
	"precision mediump float;\n"
	"uniform mat4 proj;\n"
	"attribute vec2 position;\n"
	"attribute vec2 texcoord;\n"
	"uniform float x_step;\n"
	"uniform float y_step;\n"
	"uniform ivec2 rear_screen;\n"
	"varying vec2 texcoord0;\n"
	"varying vec2 texcoord1;\n"
	"varying vec2 texcoord2;\n"
	"varying vec2 texcoord3;\n"
	"void main() {\n"
	"  float tx0 = texcoord.x + x_step;\n"
	"  float tx1 = texcoord.x + x_step * 0.36;\n"
	"  float tx2 = texcoord.x - x_step * 0.36;\n"
	"  float tx3 = texcoord.x - x_step;\n"
	"  if(rear_screen.x == 1 && rear_screen.y == 0) {\n" //left
	"    tx0 = min(0.50, tx0);\n"
	"    tx1 = min(0.50, tx1);\n"
	"  } else if(rear_screen.x == 1 && rear_screen.y == 1) {\n" //right
	"    tx2 = max(0.50, tx2);\n"
	"    tx3 = max(0.50, tx3);\n"
	"  }\n"
	"  texcoord0 = vec2(tx0, texcoord.y + y_step * 0.36);\n"
	"  texcoord1 = vec2(tx1, texcoord.y - y_step);\n"
	"  texcoord2 = vec2(tx2, texcoord.y + y_step);\n"
	"  texcoord3 = vec2(tx3, texcoord.y - y_step * 0.36);\n"
	"  gl_Position = proj * vec4(position, 0.0, 1.0);\n"
	"}\n";
static const char texture_fragment_shader_internal_blur[] =
	"precision mediump float;\n"
	"uniform sampler2D source_texture;\n"
	"varying vec2 texcoord0;\n"
	"varying vec2 texcoord1;\n"
	"varying vec2 texcoord2;\n"
	"varying vec2 texcoord3;\n"
	"void main() {\n"
	"vec4 sourceColor = ( texture2D(source_texture, texcoord0) + \n"
	"                     texture2D(source_texture, texcoord1) + "
	"                     texture2D(source_texture, texcoord2) + "
	"                     texture2D(source_texture, texcoord3) ) * 0.25;\n"
	"gl_FragColor = sourceColor;\n"
	;
static const char texture_vertex_shader_fastblur[] =
	"uniform mat4 proj;\n"
	"attribute vec2 position;\n"
	"attribute vec2 texcoord;\n"
	"varying vec2 v_texcoord;\n"
	"void main()\n"
	"{\n"
	"   gl_Position = proj * vec4(position, 0.0, 1.0);\n"
	"   v_texcoord = texcoord;\n"
	"}\n";
static const char texture_fragment_shader_fastblur[] =
	"precision mediump float;\n"
	"uniform sampler2D source0;\n"
	"uniform sampler2D source1;\n"
	"uniform sampler2D source2;\n"
	"uniform sampler2D source3;\n"
	"uniform sampler2D source4;\n"
	"uniform sampler2D source5;\n"
	"uniform float weight1;\n"
	"uniform float weight2;\n"
	"uniform float weight3;\n"
	"uniform float weight4;\n"
	"uniform float weight5;\n"
	"varying vec2 v_texcoord;\n"
	"uniform vec4 blur_region;\n"
	"uniform float alpha;\n"
	"uniform vec3 gamma;\n"
	"uniform vec3 shade;\n"
	"bool in_blur_region()\n"
	"{\n"
	"  if(v_texcoord.x >= blur_region.x && v_texcoord.x <= blur_region.y && \n"
	"      v_texcoord.y >= blur_region.z && v_texcoord.y <= blur_region.w) {\n"
	"    return true;\n"
	"  } else { \n"
	"    return false;\n"
	"  }\n"
	"}\n"
	"void main() {\n"
	  "vec4 col;\n"
	"if (in_blur_region()) {"
	  "vec4 blurCol;\n"
	"  col += texture2D(source1, v_texcoord) * weight1;\n"
	"  col += texture2D(source2, v_texcoord) * weight2;\n"
	"  col += texture2D(source3, v_texcoord) * weight3;\n"
	"  col += texture2D(source4, v_texcoord) * weight4;\n"
	"  col += texture2D(source5, v_texcoord) * weight5;\n"
	"} else {\n"
	  "col = texture2D(source0, v_texcoord);\n"
	"}\n"
	"col.rgb *= shade;\n"
	"col.rgb = pow(col.rgb, gamma);\n"
	"gl_FragColor = alpha * col;\n"
	;
#endif
#if QT_GAUSSIAN_BLUR
static const char texture_vertex_shader_gaussianblur[] =
	"uniform mat4 proj;\n"
	"attribute vec2 position;\n"
	"attribute vec2 texcoord;\n"
	"varying vec2 v_texcoord;\n"
	"uniform int samples;\n"
	"uniform vec2 tC_uni;\n"
	"uniform vec2 tL_uni[15];\n"
	"uniform vec2 tR_uni[15];\n"
	"uniform float spread;\n"
	"uniform vec2 dirstep;\n"
	"varying vec2 tC;\n"
	"varying vec2 tL[15];\n"
	"varying vec2 tR[15];\n"
	"void main()\n"
	"{\n"
	"   gl_Position = proj * vec4(position, 0.0, 1.0);\n"
	"   v_texcoord = texcoord;\n"
	"   if(mod(float(samples),2.0) == 1.0) {\n"
	"      tC = texcoord;\n"
	"      if(tC_uni.y > 0.0)\n"
	"         tC += spread * dirstep * float(tC_uni.x);\n"
	"   }\n"
	"   for(int i=0; i<samples/2; ++i) {\n"
	"      tL[i] = texcoord;\n"
	"      if(tL_uni[i].y > 0.0)\n"
	"         tL[i] += spread * dirstep * float(tL_uni[i].x);\n"
	"      tR[i] = texcoord;\n"
	"      if(tR_uni[i].x > 0.0)\n"
	"         tR[i] += spread * dirstep * float(tR_uni[i].x);\n"
	"   }\n"
	"}\n";
static const char texture_fragment_shader_gaussianblur[] =
	"precision mediump float;\n"
	"uniform sampler2D source_texture;\n"
	"varying vec2 v_texcoord;\n"
	"uniform vec4 blur_region;\n"
	"uniform float alpha;\n"
	"uniform vec3 gamma;\n"
	"uniform vec3 shade;\n"
	"uniform int samples;\n"
	"uniform vec2 tC_uni;\n"
	"uniform vec2 tL_uni[15];\n"
	"uniform vec2 tR_uni[15];\n"
	"varying vec2 tC;\n"
	"varying vec2 tL[15];\n"
	"varying vec2 tR[15];\n"
	"bool in_blur_region()\n"
	"{\n"
	"  if(v_texcoord.x >= blur_region.x && v_texcoord.x <= blur_region.y && \n"
	"      v_texcoord.y >= blur_region.z && v_texcoord.y <= blur_region.w) {\n"
	"    return true;\n"
	"    } else { \n"
	"        return false;\n"
	"    }\n"
	"}\n"
	"void main()\n"
	"{\n"
	"    vec4 col;\n"
	"    float sum = 0.0;\n"
	"    if (in_blur_region()) {"
	"        for(int i=0; i<samples/2; ++i) {\n"
	"          sum += tL_uni[i].y;\n"
	"          sum += tR_uni[i].y;\n"
	"        }\n"
	"        if(mod(float(samples),2.0) == 1.0) {\n"
	"          sum += tC_uni.y;\n"
	"          col += float(tC_uni.y/sum) * texture2D(source_texture, tC);\n"
	"        }\n"
	"        for(int i=0; i<samples/2; ++i) {\n"
	"          col += float(tL_uni[i].y/sum) * texture2D(source_texture, tL[i]);\n"
	"          col += float(tR_uni[i].y/sum) * texture2D(source_texture, tR[i]);\n"
	"        }\n"
	"    } else {\n"
	        "col = texture2D(source_texture, v_texcoord);\n"
	"    }\n"
	"    col.rgb *= shade;\n"
	"    col.rgb = pow(col.rgb, gamma);\n"
	"    gl_FragColor = alpha * col;\n"
	;
#endif
#else
static const char vertex_shader[] =
	"uniform mat4 proj;\n"
	"attribute vec2 position;\n"
	"attribute vec2 texcoord;\n"
	"varying vec2 v_texcoord;\n"
	"void main()\n"
	"{\n"
	"   gl_Position = proj * vec4(position, 0.0, 1.0);\n"
	"   v_texcoord = texcoord;\n"
	"}\n";
#endif // HKMC_CCIC
/* Declare common fragment shader uniforms */
#define FRAGMENT_CONVERT_YUV						\
	"  y *= alpha;\n"						\
	"  u *= alpha;\n"						\
	"  v *= alpha;\n"						\
	"  gl_FragColor.r = y + 1.59602678 * v;\n"			\
	"  gl_FragColor.g = y - 0.39176229 * u - 0.81296764 * v;\n"	\
	"  gl_FragColor.b = y + 2.01723214 * u;\n"			\
	"  gl_FragColor.a = alpha;\n"
static const char fragment_debug[] =
	"  gl_FragColor = vec4(0.0, 0.3, 0.0, 0.2) + gl_FragColor * 0.8;\n";
static const char fragment_brace[] =
	"}\n";
static const char texture_fragment_shader_rgba[] =
	"precision mediump float;\n"
	"varying vec2 v_texcoord;\n"
	"uniform sampler2D tex;\n"
	"uniform float alpha;\n"
	"void main()\n"
	"{\n"
	"   gl_FragColor = alpha * texture2D(tex, v_texcoord)\n;"
	;
static const char texture_fragment_shader_rgbx[] =
	"precision mediump float;\n"
	"varying vec2 v_texcoord;\n"
	"uniform sampler2D tex;\n"
	"uniform float alpha;\n"
	"void main()\n"
	"{\n"
	"   gl_FragColor.rgb = alpha * texture2D(tex, v_texcoord).rgb\n;"
	"   gl_FragColor.a = alpha;\n"
	;
static const char texture_fragment_shader_egl_external[] =
	"#extension GL_OES_EGL_image_external : require\n"
	"precision mediump float;\n"
	"varying vec2 v_texcoord;\n"
	"uniform samplerExternalOES tex;\n"
	"uniform float alpha;\n"
#if 1 // nvidia
	"uniform vec3 gamma;\n"
	"uniform vec3 shade;\n"
#endif
	"void main()\n"
	"{\n"
#if 1 // nvidia
	"   vec4 col = texture2D(tex, v_texcoord);\n"
	"   col.rgb *= shade;\n"
	"   col.rgb = pow(col.rgb, gamma);\n"
	"   gl_FragColor = alpha * col;\n"
#else
    "   gl_FragColor = alpha * texture2D(tex, v_texcoord)\n;"
#endif
	;
static const char texture_fragment_shader_y_uv[] =
	"precision mediump float;\n"
	"uniform sampler2D tex;\n"
	"uniform sampler2D tex1;\n"
	"varying vec2 v_texcoord;\n"
	"uniform float alpha;\n"
	"void main() {\n"
	"  float y = 1.16438356 * (texture2D(tex, v_texcoord).x - 0.0625);\n"
	"  float u = texture2D(tex1, v_texcoord).r - 0.5;\n"
	"  float v = texture2D(tex1, v_texcoord).g - 0.5;\n"
	FRAGMENT_CONVERT_YUV
	;
static const char texture_fragment_shader_y_u_v[] =
	"precision mediump float;\n"
	"uniform sampler2D tex;\n"
	"uniform sampler2D tex1;\n"
	"uniform sampler2D tex2;\n"
	"varying vec2 v_texcoord;\n"
	"uniform float alpha;\n"
	"void main() {\n"
	"  float y = 1.16438356 * (texture2D(tex, v_texcoord).x - 0.0625);\n"
	"  float u = texture2D(tex1, v_texcoord).x - 0.5;\n"
	"  float v = texture2D(tex2, v_texcoord).x - 0.5;\n"
	FRAGMENT_CONVERT_YUV
	;
static const char texture_fragment_shader_y_xuxv[] =
	"precision mediump float;\n"
	"uniform sampler2D tex;\n"
	"uniform sampler2D tex1;\n"
	"varying vec2 v_texcoord;\n"
	"uniform float alpha;\n"
	"void main() {\n"
	"  float y = 1.16438356 * (texture2D(tex, v_texcoord).x - 0.0625);\n"
	"  float u = texture2D(tex1, v_texcoord).g - 0.5;\n"
	"  float v = texture2D(tex1, v_texcoord).a - 0.5;\n"
	FRAGMENT_CONVERT_YUV
	;
static const char texture_fragment_shader_xyuv[] =
	"precision mediump float;\n"
	"uniform sampler2D tex;\n"
	"varying vec2 v_texcoord;\n"
	"uniform float alpha;\n"
	"void main() {\n"
	"  float y = 1.16438356 * (texture2D(tex, v_texcoord).b - 0.0625);\n"
	"  float u = texture2D(tex, v_texcoord).g - 0.5;\n"
	"  float v = texture2D(tex, v_texcoord).r - 0.5;\n"
	FRAGMENT_CONVERT_YUV
	;
static const char solid_fragment_shader[] =
	"precision mediump float;\n"
	"uniform vec4 color;\n"
	"uniform float alpha;\n"
	"void main()\n"
	"{\n"
	"   gl_FragColor = alpha * color\n;"
	;
static int
compile_shader(GLenum type, int count, const char **sources)
{
	GLuint s;
	char msg[512];
	GLint status;
	s = glCreateShader(type);
	glShaderSource(s, count, sources, NULL);
	glCompileShader(s);
	glGetShaderiv(s, GL_COMPILE_STATUS, &status);
	if (!status) {
		glGetShaderInfoLog(s, sizeof msg, NULL, msg);
		weston_log("shader info: %s\n", msg);
		return GL_NONE;
	}
	return s;
}
static int
shader_init(struct gl_shader *shader, struct gl_renderer *renderer,
		   const char *vertex_source, const char *fragment_source)
{
	char msg[512];
	GLint status;
	int count;
	const char *sources[3];
	shader->vertex_shader =
		compile_shader(GL_VERTEX_SHADER, 1, &vertex_source);
	if (renderer->fragment_shader_debug) {
		sources[0] = fragment_source;
		sources[1] = fragment_debug;
		sources[2] = fragment_brace;
		count = 3;
	} else {
		sources[0] = fragment_source;
		sources[1] = fragment_brace;
		count = 2;
	}
	shader->fragment_shader =
		compile_shader(GL_FRAGMENT_SHADER, count, sources);
	shader->program = glCreateProgram();
	glAttachShader(shader->program, shader->vertex_shader);
	glAttachShader(shader->program, shader->fragment_shader);
	glBindAttribLocation(shader->program, 0, "position");
	glBindAttribLocation(shader->program, 1, "texcoord");
	glLinkProgram(shader->program);
	glGetProgramiv(shader->program, GL_LINK_STATUS, &status);
	if (!status) {
		glGetProgramInfoLog(shader->program, sizeof msg, NULL, msg);
		weston_log("link info: %s\n", msg);
		return -1;
	}
	shader->proj_uniform = glGetUniformLocation(shader->program, "proj");
	shader->tex_uniforms[0] = glGetUniformLocation(shader->program, "tex");
	shader->tex_uniforms[1] = glGetUniformLocation(shader->program, "tex1");
	shader->tex_uniforms[2] = glGetUniformLocation(shader->program, "tex2");
	shader->alpha_uniform = glGetUniformLocation(shader->program, "alpha");
	shader->color_uniform = glGetUniformLocation(shader->program, "color");
#if 1 // nvidia
	shader->gamma_uniform = glGetUniformLocation(shader->program, "gamma");
	shader->shade_uniform = glGetUniformLocation(shader->program, "shade");
#endif
#if HKMC_CCIC
#if TILT_ENABLE
	shader->tilt_uniform = glGetUniformLocation(shader->program, "tilt");
#endif
	shader->output_size_uniform = glGetUniformLocation(shader->program, "output_size");
#if QT_FAST_BLUR
	shader->x_step_uniform = glGetUniformLocation(shader->program, "x_step");
	shader->y_step_uniform = glGetUniformLocation(shader->program, "y_step");
	shader->blur_region_uniform = glGetUniformLocation(shader->program, "blur_region");
	shader->source_texture_uniform = glGetUniformLocation(shader->program, "source_texture");
	shader->source_uniform[0] = glGetUniformLocation(shader->program, "source0");
	shader->source_uniform[1] = glGetUniformLocation(shader->program, "source1");
	shader->source_uniform[2] = glGetUniformLocation(shader->program, "source2");
	shader->source_uniform[3] = glGetUniformLocation(shader->program, "source3");
	shader->source_uniform[4] = glGetUniformLocation(shader->program, "source4");
	shader->source_uniform[5] = glGetUniformLocation(shader->program, "source5");
	shader->weight_uniform[0] = glGetUniformLocation(shader->program, "weight1");
	shader->weight_uniform[1] = glGetUniformLocation(shader->program, "weight2");
	shader->weight_uniform[2] = glGetUniformLocation(shader->program, "weight3");
	shader->weight_uniform[3] = glGetUniformLocation(shader->program, "weight4");
	shader->weight_uniform[4] = glGetUniformLocation(shader->program, "weight5");
	shader->rear_screen_uniform = glGetUniformLocation(shader->program, "rear_screen");
#elif QT_GAUSSIAN_BLUR
	shader->blur_region_uniform = glGetUniformLocation(shader->program, "blur_region");
	shader->source_texture_uniform = glGetUniformLocation(shader->program, "source_texture");
	shader->spread_uniform = glGetUniformLocation(shader->program, "spread");
	shader->dirstep_uniform = glGetUniformLocation(shader->program, "dirstep");
	shader->samples_uniform = glGetUniformLocation(shader->program, "samples");
	shader->tC_uniform = glGetUniformLocation(shader->program, "tC_uni");
	int i;
	for(i = 0; i < GAUSS_ARRAY_LENGTH; ++i) {
		char tL[11];
		char tR[11];
		snprintf(tL, sizeof tL, "tL_uni[%d]", i);
		snprintf(tR, sizeof tR, "tR_uni[%d]", i);
		shader->tL_uniform[i] = glGetUniformLocation(shader->program, tL);
		shader->tR_uniform[i] = glGetUniformLocation(shader->program, tR);
	}
	shader->rear_screen_uniform = glGetUniformLocation(shader->program, "rear_screen");
#endif
#endif // HKMC_CCIC
	return 0;
}
static void
shader_release(struct gl_shader *shader)
{
	glDeleteShader(shader->vertex_shader);
	glDeleteShader(shader->fragment_shader);
	glDeleteProgram(shader->program);
	shader->vertex_shader = 0;
	shader->fragment_shader = 0;
	shader->program = 0;
}
static void
log_extensions(const char *name, const char *extensions)
{
	const char *p, *end;
	int l;
	int len;
	l = weston_log("%s:", name);
	p = extensions;
	while (*p) {
		end = strchrnul(p, ' ');
		len = end - p;
		if (l + len > 78)
			l = weston_log_continue("\n" STAMP_SPACE "%.*s",
						len, p);
		else
			l += weston_log_continue(" %.*s", len, p);
		for (p = end; isspace(*p); p++)
			;
	}
	weston_log_continue("\n");
}
static void
log_egl_info(EGLDisplay egldpy)
{
	const char *str;
	str = eglQueryString(egldpy, EGL_VERSION);
	weston_log("EGL version: %s\n", str ? str : "(null)");
	str = eglQueryString(egldpy, EGL_VENDOR);
	weston_log("EGL vendor: %s\n", str ? str : "(null)");
	str = eglQueryString(egldpy, EGL_CLIENT_APIS);
	weston_log("EGL client APIs: %s\n", str ? str : "(null)");
	str = eglQueryString(egldpy, EGL_EXTENSIONS);
	log_extensions("EGL extensions", str ? str : "(null)");
}
static void
log_gl_info(void)
{
	const char *str;
	str = (char *)glGetString(GL_VERSION);
	weston_log("GL version: %s\n", str ? str : "(null)");
	str = (char *)glGetString(GL_SHADING_LANGUAGE_VERSION);
	weston_log("GLSL version: %s\n", str ? str : "(null)");
	str = (char *)glGetString(GL_VENDOR);
	weston_log("GL vendor: %s\n", str ? str : "(null)");
	str = (char *)glGetString(GL_RENDERER);
	weston_log("GL renderer: %s\n", str ? str : "(null)");
	str = (char *)glGetString(GL_EXTENSIONS);
	log_extensions("GL extensions", str ? str : "(null)");
}
static void
gl_renderer_output_set_border(struct weston_output *output,
			      enum gl_renderer_border_side side,
			      int32_t width, int32_t height,
			      int32_t tex_width, unsigned char *data)
{
	struct gl_output_state *go = get_output_state(output);
	if (go->borders[side].width != width ||
	    go->borders[side].height != height)
		/* In this case, we have to blow everything and do a full
		 * repaint. */
		go->border_status |= BORDER_SIZE_CHANGED | BORDER_ALL_DIRTY;
	if (data == NULL) {
		width = 0;
		height = 0;
	}
	go->borders[side].width = width;
	go->borders[side].height = height;
	go->borders[side].tex_width = tex_width;
	go->borders[side].data = data;
	go->border_status |= 1 << side;
}
static int
gl_renderer_setup(struct weston_compositor *ec, EGLSurface egl_surface);
static EGLSurface
gl_renderer_create_window_surface(struct gl_renderer *gr,
				  EGLNativeWindowType window_for_legacy,
				  void *window_for_platform,
				  const uint32_t *drm_formats,
				  unsigned drm_formats_count)
{
	EGLSurface egl_surface = EGL_NO_SURFACE;
	EGLConfig egl_config;
	egl_config = gl_renderer_get_egl_config(gr, EGL_WINDOW_BIT,
						drm_formats, drm_formats_count);
	if (egl_config == EGL_NO_CONFIG_KHR)
		return EGL_NO_SURFACE;
	log_egl_config_info(gr->egl_display, egl_config);
	if (gr->create_platform_window)
		egl_surface = gr->create_platform_window(gr->egl_display,
							 egl_config,
							 window_for_platform,
							 NULL);
	else
		egl_surface = eglCreateWindowSurface(gr->egl_display,
						     egl_config,
						     window_for_legacy, NULL);
	return egl_surface;
}
#if 1 // nvidia
static EGLSurface
gl_renderer_create_stream_surface(struct gl_renderer *gr,
				  uint32_t plane_id,
				  uint32_t crtc_id,
				  EGLint width, EGLint height,
				  EGLStreamKHR *egl_stream)
{
	EGLint stream_attribs[] = {
		EGL_STREAM_FIFO_LENGTH_KHR, 0, /* Mailbox mode */
#ifdef EGL_EXT_stream_acquire_mode
		EGL_CONSUMER_AUTO_ACQUIRE_EXT, EGL_FALSE,
#endif
		EGL_NONE
	};
	EGLAttrib output_attribs[3];
	EGLint stream_producer_attribs[] = {
		EGL_WIDTH, width,
		EGL_HEIGHT, height,
		EGL_NONE
	};
	EGLint num_layers;
	EGLOutputLayerEXT output_layer;
	EGLSurface egl_surface = EGL_NO_SURFACE;
	*egl_stream = gr->create_stream(gr->egl_display, stream_attribs);
	if (*egl_stream == EGL_NO_STREAM_KHR) {
		weston_log("Failed to create EGL stream.\n");
		goto err_egl_create_surf_base;
	}
	if (plane_id != ~0u) {
		output_attribs[0] = EGL_DRM_PLANE_EXT;
		output_attribs[1] = plane_id;
	} else {
		assert(crtc_id != ~0u);
		output_attribs[0] = EGL_DRM_CRTC_EXT;
		output_attribs[1] = crtc_id;
	}
	output_attribs[2] = EGL_NONE;
	if (gr->get_output_layers(gr->egl_display,
				  output_attribs,
				  &output_layer,
				  1, &num_layers) != EGL_TRUE) {
		weston_log("Failed to get output layer.\n");
		goto err_egl_create_surf_stream;
	}
	if (num_layers < 1) {
		weston_log("Unable to find output layers.\n");
		goto err_egl_create_surf_stream;
	}
	if (gr->stream_consumer_output(gr->egl_display, *egl_stream,
				       output_layer) != EGL_TRUE) {
		weston_log("Failed to set EGL stream consumer.\n");
		goto err_egl_create_surf_stream;
	}
	egl_surface = gr->create_stream_producer_surface(
						gr->egl_display,
						gr->egl_config,
						*egl_stream,
						stream_producer_attribs);
	if (egl_surface == EGL_NO_SURFACE) {
		weston_log("Failed to create EGL producer surface.\n");
		goto err_egl_create_surf_stream;
	}
	return egl_surface;
err_egl_create_surf_stream:
	gr->destroy_stream(gr->egl_display, *egl_stream);
	*egl_stream = EGL_NO_STREAM_KHR;
err_egl_create_surf_base:
	gl_renderer_print_egl_error_state();
	return EGL_NO_SURFACE;
}
#endif // nvidia
static int
gl_renderer_output_create(struct weston_output *output,
			  EGLSurface surface,
			  EGLStreamKHR stream)
{
	struct gl_output_state *go;
	int i;
	go = zalloc(sizeof *go);
	if (go == NULL)
		return -1;
	go->egl_surface = surface;
#if 1 // nvidia
	go->egl_stream = stream;
#endif
	for (i = 0; i < BUFFER_DAMAGE_COUNT; i++)
		pixman_region32_init(&go->buffer_damage[i]);
	wl_list_init(&go->timeline_render_point_list);
	go->begin_render_sync = EGL_NO_SYNC_KHR;
	go->end_render_sync = EGL_NO_SYNC_KHR;
	output->renderer_state = go;
#if 1 // nvidia
	go->gamma[0] = 1.0;
	go->gamma[1] = 1.0;
	go->gamma[2] = 1.0;
	go->shade[0] = 1.0;
	go->shade[1] = 1.0;
	go->shade[2] = 1.0;
#endif
	return 0;
}
static int
gl_renderer_output_window_create(struct weston_output *output,
				 EGLNativeWindowType window_for_legacy,
				 void *window_for_platform,
				 const uint32_t *drm_formats,
				 unsigned drm_formats_count)
{
	struct weston_compositor *ec = output->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	EGLSurface egl_surface = EGL_NO_SURFACE;
	int ret = 0;
	egl_surface = gl_renderer_create_window_surface(gr,
							window_for_legacy,
							window_for_platform,
							drm_formats,
							drm_formats_count);
	if (egl_surface == EGL_NO_SURFACE) {
		weston_log("failed to create egl surface\n");
		return -1;
	}
#if 1
	ret = gl_renderer_output_create(output, egl_surface, EGL_NO_STREAM_KHR);
#else
    ret = gl_renderer_output_create(output, egl_surface);
#endif
	if (ret < 0)
		weston_platform_destroy_egl_surface(gr->egl_display, egl_surface);
	return ret;
}
static int
gl_renderer_output_pbuffer_create(struct weston_output *output,
				  int width,
				  int height,
				  const uint32_t *drm_formats,
				  unsigned drm_formats_count)
{
	struct gl_renderer *gr = get_renderer(output->compositor);
	EGLConfig pbuffer_config;
	EGLSurface egl_surface;
	int ret;
	EGLint pbuffer_attribs[] = {
		EGL_WIDTH, width,
		EGL_HEIGHT, height,
		EGL_NONE
	};
	pbuffer_config = gl_renderer_get_egl_config(gr, EGL_PBUFFER_BIT,
						    drm_formats,
						    drm_formats_count);
	if (pbuffer_config == EGL_NO_CONFIG_KHR) {
		weston_log("failed to choose EGL config for PbufferSurface\n");
		return -1;
	}
	log_egl_config_info(gr->egl_display, pbuffer_config);
	egl_surface = eglCreatePbufferSurface(gr->egl_display, pbuffer_config,
					      pbuffer_attribs);
	if (egl_surface == EGL_NO_SURFACE) {
		weston_log("failed to create egl surface\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
#if 1 // nvidia
	ret = gl_renderer_output_create(output, egl_surface, EGL_NO_STREAM_KHR);
#else
    ret = gl_renderer_output_create(output, egl_surface);
#endif
	if (ret < 0)
		eglDestroySurface(gr->egl_display, egl_surface);
	return ret;
}
#if 1 // nvidia
static int
gl_renderer_output_stream_create(struct weston_output *output,
				 uint32_t plane_id, uint32_t crtc_id)
{
	struct weston_compositor *ec = output->compositor;
	struct gl_renderer *gr = get_renderer(ec);
	EGLSurface egl_surface = EGL_NO_SURFACE;
	EGLStreamKHR egl_stream = EGL_NO_STREAM_KHR;
	int ret;
	egl_surface =
		gl_renderer_create_stream_surface(gr,
						  plane_id, crtc_id,
						  output->current_mode->width,
						  output->current_mode->height,
						  &egl_stream);
	if (egl_surface == EGL_NO_SURFACE)
		return -1;
	ret = gl_renderer_output_create(output, egl_surface, egl_stream);
	if (ret < 0) {
		eglDestroySurface(gr->egl_display, egl_surface);
		gr->destroy_stream(gr->egl_display, egl_stream);
	}
	return ret;
}
#endif
static void
gl_renderer_output_destroy(struct weston_output *output)
{
	struct gl_renderer *gr = get_renderer(output->compositor);
	struct gl_output_state *go = get_output_state(output);
	struct timeline_render_point *trp, *tmp;
	int i;
	for (i = 0; i < 2; i++)
		pixman_region32_fini(&go->buffer_damage[i]);
	eglMakeCurrent(gr->egl_display,
		       EGL_NO_SURFACE, EGL_NO_SURFACE,
		       EGL_NO_CONTEXT);
	weston_platform_destroy_egl_surface(gr->egl_display, go->egl_surface);
	if (!wl_list_empty(&go->timeline_render_point_list))
		weston_log("warning: discarding pending timeline render"
			   "objects at output destruction");
	wl_list_for_each_safe(trp, tmp, &go->timeline_render_point_list, link)
		timeline_render_point_destroy(trp);
	if (go->begin_render_sync != EGL_NO_SYNC_KHR)
		gr->destroy_sync(gr->egl_display, go->begin_render_sync);
	if (go->end_render_sync != EGL_NO_SYNC_KHR)
		gr->destroy_sync(gr->egl_display, go->end_render_sync);
#if 1 // nvidia
	if (go->egl_stream != EGL_NO_STREAM_KHR)
		gr->destroy_stream(gr->egl_display, go->egl_stream);
#endif
#if HKMC_CCIC
	if (output->enable_fbo) {
#if QT_FAST_BLUR
		for (i = 0; i < FBO_REPAINT_MAX; i++) {
			glDeleteTextures(1, &go->blur_tex[i]);
			glDeleteFramebuffers(1, &go->blur_fbo[i]);
		}
#elif QT_GAUSSIAN_BLUR
		for (i = 0; i < FBO_REPAINT_MAX; i++) {
			glDeleteTextures(1, &go->blur_tex[i]);
			glDeleteFramebuffers(1, &go->blur_fbo[i]);
		}
#endif
		output->enable_fbo = false;
    }
#endif
	free(go);
}
static int
gl_renderer_create_fence_fd(struct weston_output *output)
{
	struct gl_output_state *go = get_output_state(output);
	struct gl_renderer *gr = get_renderer(output->compositor);
	int fd;
	if (go->end_render_sync == EGL_NO_SYNC_KHR)
		return -1;
	fd = gr->dup_native_fence_fd(gr->egl_display, go->end_render_sync);
	if (fd == EGL_NO_NATIVE_FENCE_FD_ANDROID)
		return -1;
	return fd;
}
static void
gl_renderer_destroy(struct weston_compositor *ec)
{
	struct gl_renderer *gr = get_renderer(ec);
	struct dmabuf_image *image, *next;
	wl_signal_emit(&gr->destroy_signal, gr);
	if (gr->has_bind_display)
		gr->unbind_display(gr->egl_display, ec->wl_display);
	/* Work around crash in egl_dri2.c's dri2_make_current() - when does this apply? */
	eglMakeCurrent(gr->egl_display,
		       EGL_NO_SURFACE, EGL_NO_SURFACE,
		       EGL_NO_CONTEXT);
	wl_list_for_each_safe(image, next, &gr->dmabuf_images, link)
		dmabuf_image_destroy(image);
	if (gr->dummy_surface != EGL_NO_SURFACE)
		weston_platform_destroy_egl_surface(gr->egl_display,
						    gr->dummy_surface);
	eglTerminate(gr->egl_display);
	eglReleaseThread();
	wl_list_remove(&gr->output_destroy_listener.link);
	wl_array_release(&gr->vertices);
	wl_array_release(&gr->vtxcnt);
	if (gr->fragment_binding)
		weston_binding_destroy(gr->fragment_binding);
	if (gr->fan_binding)
		weston_binding_destroy(gr->fan_binding);
	free(gr);
}
/** Checks whether a platform EGL client extension is supported
 *
 * \param ec The weston compositor
 * \param extension_suffix The EGL client extension suffix
 * \return 1 if supported, 0 if using fallbacks, -1 unsupported
 *
 * This function checks whether a specific platform_* extension is supported
 * by EGL.
 *
 * The extension suffix should be the suffix of the platform extension (that
 * specifies a platform argument as defined in EGL_EXT_platform_base). For
 * example, passing "foo" will check whether either "EGL_KHR_platform_foo",
 * "EGL_EXT_platform_foo", or "EGL_MESA_platform_foo" is supported.
 *
 * The return value is 1:
 *   - if the supplied EGL client extension is supported.
 * The return value is 0:
 *   - if the platform_base client extension isn't supported so will
 *     fallback to eglGetDisplay and friends.
 * The return value is -1:
 *   - if the supplied EGL client extension is not supported.
 */
static int
gl_renderer_supports(struct weston_compositor *ec,
		     const char *extension_suffix)
{
	static const char *extensions = NULL;
	char s[64];
	if (!extensions) {
		extensions = (const char *) eglQueryString(
			EGL_NO_DISPLAY, EGL_EXTENSIONS);
		if (!extensions)
			return 0;
		log_extensions("EGL client extensions",
			       extensions);
	}
	if (!weston_check_egl_extension(extensions, "EGL_EXT_platform_base"))
		return 0;
	snprintf(s, sizeof s, "EGL_KHR_platform_%s", extension_suffix);
	if (weston_check_egl_extension(extensions, s))
		return 1;
	snprintf(s, sizeof s, "EGL_EXT_platform_%s", extension_suffix);
	if (weston_check_egl_extension(extensions, s))
		return 1;
	snprintf(s, sizeof s, "EGL_MESA_platform_%s", extension_suffix);
	if (weston_check_egl_extension(extensions, s))
		return 1;
	/* at this point we definitely have some platform extensions but
	 * haven't found the supplied platform, so chances are it's
	 * not supported. */
	return -1;
}
static const char *
platform_to_extension(EGLenum platform)
{
	switch (platform) {
	case EGL_PLATFORM_GBM_KHR:
		return "gbm";
	case EGL_PLATFORM_WAYLAND_KHR:
		return "wayland";
	case EGL_PLATFORM_X11_KHR:
		return "x11";
	case EGL_PLATFORM_SURFACELESS_MESA:
		return "surfaceless";
#if 1 // nvidia
	case EGL_PLATFORM_DEVICE_EXT:
		return "device";
#endif
	default:
		assert(0 && "bad EGL platform enum");
	}
}
static void
output_handle_destroy(struct wl_listener *listener, void *data)
{
	struct gl_renderer *gr;
	struct weston_output *output = data;
	gr = container_of(listener, struct gl_renderer,
			  output_destroy_listener);
	if (wl_list_empty(&output->compositor->output_list))
		eglMakeCurrent(gr->egl_display, gr->dummy_surface,
			       gr->dummy_surface, gr->egl_context);
}
static int
gl_renderer_create_pbuffer_surface(struct gl_renderer *gr) {
	EGLConfig pbuffer_config;
	static const EGLint pbuffer_attribs[] = {
		EGL_WIDTH, 10,
		EGL_HEIGHT, 10,
		EGL_NONE
	};
	pbuffer_config = gl_renderer_get_egl_config(gr, EGL_PBUFFER_BIT,
						    NULL, 0);
	if (pbuffer_config == EGL_NO_CONFIG_KHR) {
		weston_log("failed to choose EGL config for PbufferSurface\n");
		return -1;
	}
	gr->dummy_surface = eglCreatePbufferSurface(gr->egl_display,
						    pbuffer_config,
						    pbuffer_attribs);
	if (gr->dummy_surface == EGL_NO_SURFACE) {
		weston_log("failed to create PbufferSurface\n");
		return -1;
	}
	return 0;
}
static int
gl_renderer_display_create(struct weston_compositor *ec,
			   EGLenum platform,
			   void *native_display,
#if 1 // nvidia
			   const EGLint *platform_attribs,
#endif
			   EGLint egl_surface_type,
			   const uint32_t *drm_formats,
			   unsigned drm_formats_count)
{
	struct gl_renderer *gr;
	EGLint major, minor;
	int supports = 0;
	if (platform) {
		supports = gl_renderer_supports(
			ec, platform_to_extension(platform));
		if (supports < 0)
			return -1;
	}
	/* Surfaceless is unusable without platform_base extension */
	if (supports == 0 && platform == EGL_PLATFORM_SURFACELESS_MESA)
		return -1;
	gr = zalloc(sizeof *gr);
	if (gr == NULL)
		return -1;
	gr->base.read_pixels = gl_renderer_read_pixels;
	gr->base.repaint_output = gl_renderer_repaint_output;
	gr->base.flush_damage = gl_renderer_flush_damage;
	gr->base.attach = gl_renderer_attach;
	gr->base.surface_set_color = gl_renderer_surface_set_color;
	gr->base.destroy = gl_renderer_destroy;
	gr->base.surface_get_content_size =
		gl_renderer_surface_get_content_size;
	gr->base.surface_copy_content = gl_renderer_surface_copy_content;
	gr->platform = platform;
#if 1 // nvidia
	gr->base.attach_eglstream_consumer =
		gl_renderer_attach_eglstream_consumer;
#endif
	gr->egl_display = NULL;
#if 1 // nvidia
	gr->base.set_output_gamma = gl_renderer_set_output_gamma;
	gr->base.set_output_shade = gl_renderer_set_output_shade;
#endif
	/* extension_suffix is supported */
	if (supports) {
		if (!get_platform_display) {
			get_platform_display = (void *) eglGetProcAddress(
					"eglGetPlatformDisplayEXT");
		}
		/* also wrap this in the supports check because
		 * eglGetProcAddress can return non-NULL and still not
		 * support the feature at runtime, so ensure the
		 * appropriate extension checks have been done. */
		if (get_platform_display && platform) {
			gr->egl_display = get_platform_display(platform,
							       native_display,
							       platform_attribs); // nvidia
		}
	}
	if (!gr->egl_display) {
		weston_log("warning: either no EGL_EXT_platform_base "
			   "support or specific platform support; "
			   "falling back to eglGetDisplay.\n");
		gr->egl_display = eglGetDisplay(native_display);
	}
	if (gr->egl_display == EGL_NO_DISPLAY) {
		weston_log("failed to create display\n");
		goto fail;
	}
	if (!eglInitialize(gr->egl_display, &major, &minor)) {
		weston_log("failed to initialize display\n");
		goto fail_with_error;
	}
	log_egl_info(gr->egl_display);
	ec->renderer = &gr->base;
	if (gl_renderer_setup_egl_extensions(ec) < 0)
		goto fail_with_error;
#if 1 // nvidia
	if (!gr->has_configless_context ||
	    platform == EGL_PLATFORM_DEVICE_EXT) {
#else
    if (!gr->has_configless_context) {
#endif
		if (!gr->has_surfaceless_context)
			egl_surface_type |= EGL_PBUFFER_BIT;
		gr->egl_config = gl_renderer_get_egl_config(gr,
							    egl_surface_type,
							    drm_formats,
							    drm_formats_count);
		if (gr->egl_config == EGL_NO_CONFIG_KHR) {
			weston_log("failed to choose EGL config\n");
			goto fail_terminate;
		}
	}
	ec->capabilities |= WESTON_CAP_ROTATION_ANY;
	ec->capabilities |= WESTON_CAP_CAPTURE_YFLIP;
	ec->capabilities |= WESTON_CAP_VIEW_CLIP_MASK;
	if (gr->has_native_fence_sync && gr->has_wait_sync)
		ec->capabilities |= WESTON_CAP_EXPLICIT_SYNC;
#if 1 // nvidia
	if (platform == EGL_PLATFORM_DEVICE_EXT) {
		if (!gr->has_egl_output_base ||
		    !gr->has_egl_output_drm ||
		    !gr->has_egl_stream ||
		    !gr->has_egl_stream_producer_eglsurface ||
		    !gr->has_egl_stream_consumer_egloutput ||
		    !gr->has_egl_stream_attrib ||
		    !gr->has_egl_stream_acquire_mode) {
			weston_log("following required extensions not "
				"supported:\n%s%s%s%s%s%s%s",
				(gr->has_egl_output_base                ? "" :
					"    EGL_EXT_output_base\n"),
				(gr->has_egl_output_drm                 ? "" :
					"    EGL_EXT_output_drm\n"),
				(gr->has_egl_stream                     ? "" :
					"    EGL_KHR_stream\n"),
				(gr->has_egl_stream_producer_eglsurface ? "" :
					"    EGL_KHR_stream_producer_eglsurface\n"),
				(gr->has_egl_stream_consumer_egloutput  ? "" :
					"    EGL_EXT_stream_consumer_egloutput\n"),
				(gr->has_egl_stream_attrib              ? "" :
					"    EGL_NV_stream_attrib\n"),
				(gr->has_egl_stream_acquire_mode        ? "" :
					"    EGL_EXT_stream_acquire_mode\n"));
			goto fail_terminate;
		}
		if (!gr->has_egl_stream_consumer_gltexture ||
		    !gr->has_egl_wayland_eglstream)
			weston_log("warning: following required extensions for "
				"EGL client frame presentation through "
				"EGLDevice not supported:\n%s%s",
				(gr->has_egl_stream_consumer_gltexture ? "" :
					"    EGL_KHR_stream_consumer_gltexture\n"),
				(gr->has_egl_wayland_eglstream         ? "" :
					"    EGL_WL_wayland_eglstream\n"));
		if (!gr->has_egl_output_drm_flip_event)
			weston_log("warning: EGL page flip event notification "
				   "not supported\n");
	}
#endif
	wl_list_init(&gr->dmabuf_images);
	if (gr->has_dmabuf_import) {
		gr->base.import_dmabuf = gl_renderer_import_dmabuf;
		gr->base.query_dmabuf_formats =
			gl_renderer_query_dmabuf_formats;
		gr->base.query_dmabuf_modifiers =
			gl_renderer_query_dmabuf_modifiers;
	}
	if (gr->has_surfaceless_context) {
		weston_log("EGL_KHR_surfaceless_context available\n");
		gr->dummy_surface = EGL_NO_SURFACE;
	} else {
		weston_log("EGL_KHR_surfaceless_context unavailable. "
			   "Trying PbufferSurface\n");
		if (gl_renderer_create_pbuffer_surface(gr) < 0)
			goto fail_with_error;
	}
	wl_display_add_shm_format(ec->wl_display, WL_SHM_FORMAT_RGB565);
	wl_display_add_shm_format(ec->wl_display, WL_SHM_FORMAT_YUV420);
	wl_display_add_shm_format(ec->wl_display, WL_SHM_FORMAT_NV12);
	wl_display_add_shm_format(ec->wl_display, WL_SHM_FORMAT_YUYV);
	wl_signal_init(&gr->destroy_signal);
	if (gl_renderer_setup(ec, gr->dummy_surface) < 0) {
		if (gr->dummy_surface != EGL_NO_SURFACE)
			weston_platform_destroy_egl_surface(gr->egl_display,
							    gr->dummy_surface);
		goto fail_with_error;
	}
	return 0;
fail_with_error:
	gl_renderer_print_egl_error_state();
fail_terminate:
	eglTerminate(gr->egl_display);
fail:
	free(gr);
	return -1;
}
static int
compile_shaders(struct weston_compositor *ec)
{
	struct gl_renderer *gr = get_renderer(ec);
#if HKMC_CCIC
#if QT_FAST_BLUR
	gr->texture_shader_internal_blur.vertex_source = texture_vertex_shader_internal_blur;
	gr->texture_shader_internal_blur.fragment_source = texture_fragment_shader_internal_blur;
	gr->texture_shader_fastblur.vertex_source = texture_vertex_shader_fastblur;
	gr->texture_shader_fastblur.fragment_source = texture_fragment_shader_fastblur;
#elif QT_GAUSSIAN_BLUR
	gr->texture_shader_gaussianblur.vertex_source = texture_vertex_shader_gaussianblur;
	gr->texture_shader_gaussianblur.fragment_source = texture_fragment_shader_gaussianblur;
#endif
#endif
	gr->texture_shader_rgba.vertex_source = vertex_shader;
	gr->texture_shader_rgba.fragment_source = texture_fragment_shader_rgba;
	gr->texture_shader_rgbx.vertex_source = vertex_shader;
	gr->texture_shader_rgbx.fragment_source = texture_fragment_shader_rgbx;
	gr->texture_shader_egl_external.vertex_source = vertex_shader;
	gr->texture_shader_egl_external.fragment_source =
		texture_fragment_shader_egl_external;
	gr->texture_shader_y_uv.vertex_source = vertex_shader;
	gr->texture_shader_y_uv.fragment_source = texture_fragment_shader_y_uv;
	gr->texture_shader_y_u_v.vertex_source = vertex_shader;
	gr->texture_shader_y_u_v.fragment_source =
		texture_fragment_shader_y_u_v;
	gr->texture_shader_y_xuxv.vertex_source = vertex_shader;
	gr->texture_shader_y_xuxv.fragment_source =
		texture_fragment_shader_y_xuxv;
	gr->texture_shader_xyuv.vertex_source = vertex_shader;
	gr->texture_shader_xyuv.fragment_source = texture_fragment_shader_xyuv;
	gr->solid_shader.vertex_source = vertex_shader;
	gr->solid_shader.fragment_source = solid_fragment_shader;
	return 0;
}
static void
fragment_debug_binding(struct weston_keyboard *keyboard,
		       const struct timespec *time,
		       uint32_t key, void *data)
{
	struct weston_compositor *ec = data;
	struct gl_renderer *gr = get_renderer(ec);
	struct weston_output *output;
	gr->fragment_shader_debug = !gr->fragment_shader_debug;
	shader_release(&gr->texture_shader_rgba);
	shader_release(&gr->texture_shader_rgbx);
	shader_release(&gr->texture_shader_egl_external);
	shader_release(&gr->texture_shader_y_uv);
	shader_release(&gr->texture_shader_y_u_v);
	shader_release(&gr->texture_shader_y_xuxv);
	shader_release(&gr->texture_shader_xyuv);
	shader_release(&gr->solid_shader);
	/* Force use_shader() to call glUseProgram(), since we need to use
	 * the recompiled version of the shader. */
	gr->current_shader = NULL;
	wl_list_for_each(output, &ec->output_list, link)
		weston_output_damage(output);
}
static void
fan_debug_repaint_binding(struct weston_keyboard *keyboard,
			  const struct timespec *time,
			  uint32_t key, void *data)
{
	struct weston_compositor *compositor = data;
	struct gl_renderer *gr = get_renderer(compositor);
	gr->fan_debug = !gr->fan_debug;
	weston_compositor_damage_all(compositor);
}
static uint32_t
get_gl_version(void)
{
	const char *version;
	int major, minor;
	version = (const char *) glGetString(GL_VERSION);
	if (version &&
	    (sscanf(version, "%d.%d", &major, &minor) == 2 ||
	     sscanf(version, "OpenGL ES %d.%d", &major, &minor) == 2)) {
		return GR_GL_VERSION(major, minor);
	}
	return GR_GL_VERSION_INVALID;
}
static int
gl_renderer_setup(struct weston_compositor *ec, EGLSurface egl_surface)
{
	struct gl_renderer *gr = get_renderer(ec);
	const char *extensions;
	EGLBoolean ret;
	EGLint context_attribs[16] = {
		EGL_CONTEXT_CLIENT_VERSION, 0,
	};
	unsigned int nattr = 2;
	if (!eglBindAPI(EGL_OPENGL_ES_API)) {
		weston_log("failed to bind EGL_OPENGL_ES_API\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	/*
	 * Being the compositor we require minimum output latency,
	 * so request a high priority context for ourselves - that should
	 * reschedule all of our rendering and its dependencies to be completed
	 * first. If the driver doesn't permit us to create a high priority
	 * context, it will fallback to the default priority (MEDIUM).
	 */
	if (gr->has_context_priority) {
		context_attribs[nattr++] = EGL_CONTEXT_PRIORITY_LEVEL_IMG;
		context_attribs[nattr++] = EGL_CONTEXT_PRIORITY_HIGH_IMG;
	}
	assert(nattr < ARRAY_LENGTH(context_attribs));
	context_attribs[nattr] = EGL_NONE;
	/* try to create an OpenGLES 3 context first */
	context_attribs[1] = 3;
	gr->egl_context = eglCreateContext(gr->egl_display, gr->egl_config,
					   EGL_NO_CONTEXT, context_attribs);
	if (gr->egl_context == NULL) {
		/* and then fallback to OpenGLES 2 */
		context_attribs[1] = 2;
		gr->egl_context = eglCreateContext(gr->egl_display,
						   gr->egl_config,
						   EGL_NO_CONTEXT,
						   context_attribs);
		if (gr->egl_context == NULL) {
			weston_log("failed to create context\n");
			gl_renderer_print_egl_error_state();
			return -1;
		}
	}
	if (gr->has_context_priority) {
		EGLint value = EGL_CONTEXT_PRIORITY_MEDIUM_IMG;
		eglQueryContext(gr->egl_display, gr->egl_context,
				EGL_CONTEXT_PRIORITY_LEVEL_IMG, &value);
		if (value != EGL_CONTEXT_PRIORITY_HIGH_IMG) {
			weston_log("Failed to obtain a high priority context.\n");
			/* Not an error, continue on as normal */
		}
	}
	ret = eglMakeCurrent(gr->egl_display, egl_surface,
			     egl_surface, gr->egl_context);
	if (ret == EGL_FALSE) {
		weston_log("Failed to make EGL context current.\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	gr->gl_version = get_gl_version();
	if (gr->gl_version == GR_GL_VERSION_INVALID) {
		weston_log("warning: failed to detect GLES version, "
			   "defaulting to 2.0.\n");
		gr->gl_version = GR_GL_VERSION(2, 0);
	}
	log_gl_info();
	gr->image_target_texture_2d =
		(void *) eglGetProcAddress("glEGLImageTargetTexture2DOES");
	extensions = (const char *) glGetString(GL_EXTENSIONS);
	if (!extensions) {
		weston_log("Retrieving GL extension string failed.\n");
		return -1;
	}
	if (!weston_check_egl_extension(extensions, "GL_EXT_texture_format_BGRA8888")) {
		weston_log("GL_EXT_texture_format_BGRA8888 not available\n");
		return -1;
	}
	if (weston_check_egl_extension(extensions, "GL_EXT_read_format_bgra"))
		ec->read_format = PIXMAN_a8r8g8b8;
	else
		ec->read_format = PIXMAN_a8b8g8r8;
	if (gr->gl_version >= GR_GL_VERSION(3, 0) ||
	    weston_check_egl_extension(extensions, "GL_EXT_unpack_subimage"))
		gr->has_unpack_subimage = true;
	if (gr->gl_version >= GR_GL_VERSION(3, 0) ||
	    weston_check_egl_extension(extensions, "GL_EXT_texture_rg"))
		gr->has_gl_texture_rg = true;
	if (weston_check_egl_extension(extensions, "GL_OES_EGL_image_external"))
		gr->has_egl_image_external = true;
	glActiveTexture(GL_TEXTURE0);
	if (compile_shaders(ec))
		return -1;
	gr->fragment_binding =
		weston_compositor_add_debug_binding(ec, KEY_S,
						    fragment_debug_binding,
						    ec);
	gr->fan_binding =
		weston_compositor_add_debug_binding(ec, KEY_F,
						    fan_debug_repaint_binding,
						    ec);
	gr->output_destroy_listener.notify = output_handle_destroy;
	wl_signal_add(&ec->output_destroyed_signal,
		      &gr->output_destroy_listener);
	weston_log("GL ES 2 renderer features:\n");
	weston_log_continue(STAMP_SPACE "read-back format: %s\n",
		ec->read_format == PIXMAN_a8r8g8b8 ? "BGRA" : "RGBA");
	weston_log_continue(STAMP_SPACE "wl_shm sub-image to texture: %s\n",
			    gr->has_unpack_subimage ? "yes" : "no");
	weston_log_continue(STAMP_SPACE "EGL Wayland extension: %s\n",
			    gr->has_bind_display ? "yes" : "no");
	return 0;
}
#if 1 // nvidia
static int
gl_renderer_get_devices(EGLint max_devices, EGLDeviceEXT *devices,
			EGLint *num_devices)
{
	const char *extensions;
	PFNEGLQUERYDEVICESEXTPROC query_devices;
	extensions = (const char *)eglQueryString(EGL_NO_DISPLAY, EGL_EXTENSIONS);
	if (!extensions) {
		weston_log("Retrieving EGL extension string failed.\n");
		return -1;
	}
	if (!weston_check_egl_extension(extensions, "EGL_EXT_device_base") &&
	    (!weston_check_egl_extension(extensions, "EGL_EXT_device_query") ||
	     !weston_check_egl_extension(extensions, "EGL_EXT_device_enumeration"))) {
		weston_log("EGL_EXT_device_base not supported\n");
		return -1;
	}
	query_devices = (void *) eglGetProcAddress("eglQueryDevicesEXT");
	if (!query_devices) {
		weston_log("Failed to get eglQueryDevicesEXT function\n");
		return -1;
	}
	if (query_devices(max_devices, devices, num_devices) != EGL_TRUE) {
		weston_log("Failed to query EGL Devices\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	return 0;
}
static int
gl_renderer_get_drm_device_file(EGLDeviceEXT device,
				const char **drm_device_file)
{
	const char *exts;
	PFNEGLQUERYDEVICESTRINGEXTPROC query_device_string;
	exts = (const char *)eglQueryString(EGL_NO_DISPLAY, EGL_EXTENSIONS);
	if (!exts) {
		weston_log("Retrieving EGL extension string failed.\n");
		return -1;
	}
	if (!weston_check_egl_extension(exts, "EGL_EXT_device_base") &&
	    (!weston_check_egl_extension(exts, "EGL_EXT_device_query") ||
	     !weston_check_egl_extension(exts, "EGL_EXT_device_enumeration"))) {
		weston_log("EGL_EXT_device_base not supported.\n");
		return -1;
	}
	query_device_string =
		(void *) eglGetProcAddress("eglQueryDeviceStringEXT");
	if (!query_device_string) {
		weston_log("Failed to get eglQueryDeviceStringEXT function\n");
		return -1;
	}
	exts = query_device_string(device, EGL_EXTENSIONS);
	if (!exts) {
		weston_log("Retrieving EGL extension string failed.\n");
		return -1;
	}
	if (!weston_check_egl_extension(exts, "EGL_EXT_device_drm")) {
		weston_log("EGL_EXT_device_drm not supported.\n");
		return -1;
	}
	(*drm_device_file) = query_device_string(device,
						 EGL_DRM_DEVICE_FILE_EXT);
	if (*drm_device_file == NULL) {
		weston_log("Failed to query DRM device name.\n");
		gl_renderer_print_egl_error_state();
		return -1;
	}
	return 0;
}
#endif // nvidia
WL_EXPORT struct gl_renderer_interface gl_renderer_interface = {
	.display_create = gl_renderer_display_create,
	.output_window_create = gl_renderer_output_window_create,
	.output_pbuffer_create = gl_renderer_output_pbuffer_create,
#if 1 // nvidia
	.output_stream_create = gl_renderer_output_stream_create,
#endif
	.output_destroy = gl_renderer_output_destroy,
	.output_set_border = gl_renderer_output_set_border,
	.create_fence_fd = gl_renderer_create_fence_fd,
#if 1 // nvidia
	.get_devices = gl_renderer_get_devices,
	.get_drm_device_file = gl_renderer_get_drm_device_file,
	.output_stream_flip = gl_renderer_output_stream_flip,
	.output_stream_enable_auto_acquire =
		gl_renderer_output_stream_enable_auto_acquire
#endif
};