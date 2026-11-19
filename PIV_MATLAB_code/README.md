# PIV analysis of organoid–collagen deformation

MATLAB code to measure the deformation of the collagen matrix around a growing
organoid from time-lapse images, using PIV (PIVlab), and to post-process the
results: Lagrangian displacements, deformation gradient and principal stretches,
interface velocities, matrix degradation velocity, and stress at the interface
for a hyperelastic collagen model. It also produces the figures and videos.

## Requirements

- **MATLAB R2023b or newer** (`smoothdata2` was introduced in R2023b; the code
  also uses `fontname`/`fontsize` (R2023a), `clim` (R2022a) and `exportgraphics`).
- **Image Processing Toolbox** (`strel`, `imopen`, `imgaussfilt`, `imbinarize`,
  `imflatfield`, `imtranslate`, `bwboundaries`, `poly2mask`, `bwdist`, ...).
- **[PIVlab](https://github.com/Shrediquette/PIVlab)** (GUI), version: `<fill in>`.

No other toolboxes are needed.

## Files

| File | Purpose |
|---|---|
| `PIV_main.m` | Main script: pre-processing, drift correction, PIV post-processing, mechanics and polar plots. Run **section by section** (see Workflow). All of its helper functions are local functions at the end of the file. |
| `field_in_band.m` | Function: band-averaged profile of a scalar/vector/tensor field along the interface (called by `PIV_main.m`). |
| `plot_rad_disp_video.m` | Video `rad_disp.mp4`: deformation of material circles around the organoid. |
| `plot_deg_videos.m` | Video `deg_mesh.mp4`: interface in the reference configuration (reference vs current), with scale bar and time stamp. |
| `plot_def_videos.m` | Video `def_mesh.mp4`: deformed PIV mesh and inclusion over the image. |
| `plot_principal_stretches_video.m` | Video `stretch_glyphs.mp4`: principal-stretch glyphs on the reference image. |

The video scripts are **scripts**, not functions: they read their inputs from the
MATLAB workspace created by `PIV_main.m`.

## Input data

Place the raw images in the working folder (the code uses relative paths):

| Pattern | Content |
|---|---|
| `z_t%03d_c002.tif` | Channel used to segment the organoid (masks) |
| `z_t%03d_c004.tif` | Collagen channel used for PIV |

`%03d` is the frame number (`001`, `002`, ...). By default there are 60 frames.

Data: `<link to the image dataset (Zenodo, etc.)>`

## Parameters

Set at the top of `PIV_main.m` (section 0.1):

| Variable | Default | Meaning |
|---|---|---|
| `offset` | 5 | Mask dilation (px) |
| `nFrames` | 60 | Number of images |
| `dx` | 0.2 | Pixel size (µm/px) |
| `dt` | 20 | Time between frames (min) |

> **Note:** `analizePIVresults` is deliberately called with
> `dt = 1`, so velocities inside `PIV_results.mat` are in px/frame. They are
> converted to µm/h later using `dx` and `dt`.

The material parameters (`params.mu`, `params.lm`, `params.m`, `params.ms`,
`params.mg`, `params.alp`) are set in the section "Compute fth and traction".

## Workflow

Run `PIV_main.m` **one section at a time** (Ctrl+Enter). It cannot be run in one
go, because two steps are done by hand in the PIVlab GUI.

1. **0.1 Initialization:** parameters.
2. **0.2 Segment organoids:** `z_t*_c002.tif` → masks `m_t*.tif`.
3. **0.3 Contrast enhancing:** `z_t*_c004.tif` → `h_t*.tif`.
4. **0.4 PIVlab, first pass (manual):**
   - Load the images `h_t*.tif` and the masks `m_t*.tif`.
   - Analyse, then apply the velocity-based validation.
   - Export all frames and save them as **`PIVlab.mat`**, containing at least
     `u_filtered` and `v_filtered`.
5. **0.5–0.6 Drift correction:** estimates the rigid drift from the field border,
   then writes drift-corrected, cropped images `dc_h_t*.tif` and masks
   `dc_m_t*.tif`.
6. **0.7 PIVlab, second pass (manual):** same as step 4, on `dc_h_t*.tif` and
   `dc_m_t*.tif`. Save it as **`dc_PIVlab.mat`**, containing `u_filtered`,
   `v_filtered`, `x` and `y`.
7. **1.0 Analyse PIV results:** `analizePIVresults` computes the Lagrangian
   trajectories `chi`, displacements `disp`, the deformation gradient `F`,
   `C = FᵀF`, principal stretches/directions `pp`, the interface in the current
   (`pn`) and reference (`p0n`) configurations, and the interface velocities
   (`vn`, `Vd`, `vCn`, `vd`). Everything is saved to **`PIV_results.mat`**.
8. **1.1 Load:** loads `PIV_results.mat` and sets `dx`, `dt`.
9. **1.2 Tight mask:** `pnPlots` (contours for the plots).
10. **2.5 Central point and normal stretch:** organoid centre `cx`, `cy` and
    radius (**required by `plot_rad_disp_video.m`**; `cx`, `cy` are not stored in
    `PIV_results.mat`).
11. **Radial displacement:** radial displacement against initial distance to the
    boundary.
12. **Test Smoothing:** deformation gradient by local least squares
    (`defgrad_pls`), giving `outSM`.
13. **Compute fth and traction:** strain energy and stress for the hyperelastic
    model, giving `fth` and `tnsn`.
14. **Polar plot:** degradation velocity and `fth` around the interface at the
    frames in `tStps`, exported as `vd_fth_pol_t%03d.pdf/.png`.


### Videos

Run after `PIV_main.m` up to at least section 2.5, in the same MATLAB session:

```matlab
plot_rad_disp_video
plot_deg_videos
plot_principal_stretches_video
plot_def_videos
```

Each script writes its `.mp4` to the working folder. The other workspace inputs
each one needs:

| Script | Needs |
|---|---|
| `plot_rad_disp_video` | `cx, cy, X, Y, chi, pn, nFrames, dt`, images `dc_h_t*.tif` |
| `plot_deg_videos` | `p0n, nFrames, dx, dt`, image `dc_h_t001.tif` |
| `plot_def_videos` | `chi, pnPlots, nFrames, dx, dt`, images `dc_h_t*.tif` |
| `plot_principal_stretches_video` | `X, Y, p0n, pp, dx, dt`, image `dc_h_t001.tif` |

The scripts create temporary variables in the workspace (for example `t`, `v`,
`I`, `F`, `H`). `F` and `H` are overwritten by `plot_principal_stretches_video`.
If you need the originals again, reload `PIV_results.mat`.

## Files generated

| File | Step |
|---|---|
| `m_t*.tif`, `h_t*.tif` | 0.2, 0.3 |
| `PIVlab.mat` | 0.4 (manual) |
| `dc_h_t*.tif`, `dc_m_t*.tif` | 0.6 |
| `dc_PIVlab.mat` | 0.7 (manual) |
| `PIV_results.mat` | 1.0 |
| `vd_fth_pol_t*.pdf`, `vd_fth_img_pol_t*.png` | Polar plot |
| `rad_disp.mp4`, `deg_mesh.mp4`, `def_mesh.mp4`, `stretch_glyphs.mp4` | Video scripts |

## Notes

- Interface velocities `vn`, `Vd`, `vd` use a Savitzky–Golay temporal window of
  `sgFrameLen = 7` frames, so they are only defined for frames 4 to `nFrames − 3`.
  In the polar plot, frames of the averaging window without data are skipped
  (with a warning).
- `C`, `pp` and `outSM` are Lagrangian: they are defined on the reference PIV
  grid `X`, `Y`.
- The variable `disp` (displacements) hides MATLAB's built-in `disp` function
  while it is in the workspace.
