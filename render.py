#!/usr/bin/env python3
"""
render.py — headless ParaView: load XDMF mesh(es), render "Surface With
Edges" in solid color, save a PNG next to each .xdmf (same base name). Uses
an isometric camera for 3-D meshes (nonzero z extent) and a top-down camera
for flat 2-D meshes.

Usage:
    pvbatch --force-offscreen-rendering render.py mesh.xdmf [width] [height]
    pvbatch --force-offscreen-rendering render.py some_folder/ [width] [height]

If a folder is given, every '*.xdmf' file directly inside it is rendered,
each to a .png sharing its own base name (e.g. meshes/kuhn.xdmf ->
meshes/kuhn.png).
"""
import sys
from pathlib import Path

from paraview.simple import *


def configure_light_kit(view) -> None:
    """
    Values below mirror the Light Inspector panel in Paraview GUI:
      Key:  Warm 0.60  Int 0.75  Ele 50   Azi 10
      Fill: Warm 0.40  K:F 3.00  Ele -75  Azi -10
      Back: Warm 0.50  K:B 1.10  Ele 0    Azi 110
      Head: Warm 0.50  K:H 3.00
      Maintain Luminance: off
    """
    view.UseLight = 1   # "Light Kit" checkbox: on

    view.KeyLightWarmth    = 0.60
    view.KeyLightIntensity = 0.75
    view.KeyLightElevation = 50
    view.KeyLightAzimuth   = 10

    view.FillLightWarmth    = 0.40
    view.FillLightKFRatio   = 3.00
    view.FillLightElevation = -75
    view.FillLightAzimuth   = -10

    view.BackLightWarmth    = 0.50
    view.BackLightKBRatio   = 1.10
    view.BackLightElevation = 0
    view.BackLightAzimuth   = 110

    view.HeadLightWarmth  = 0.50
    view.HeadLightKHRatio = 3.00

    view.MaintainLuminance = 0

    # Smooths jagged edges so mesh lines read as solid black rather than
    # dithered/grayish -- improves perceived edge-vs-fill contrast.
    view.UseFXAA = 1


def render_one(xdmf_path: Path, view, width: int, height: int) -> None:
    png_path = xdmf_path.with_suffix(".png")

    reader = OpenDataFile(str(xdmf_path))
    reader.UpdatePipeline()

    disp = Show(reader, view)
    disp.Representation = 'Surface With Edges'
    ColorBy(disp, None)   # solid color, not colored by an array
    disp.AmbientColor = [0.8, 0.8, 0.85]
    disp.DiffuseColor = [0.8, 0.8, 0.85]

    # Ambient/Diffuse balance: with only directional Light Kit lighting,
    # Raising Ambient gives every face a brightness floor;
    # lowering Diffuse slightly keeps some shading for depth without the
    # harsh face-to-face contrast.
    disp.Ambient  = 0.5
    disp.Diffuse  = 0.65
    disp.Specular = 0.0

    disp.EdgeColor = [0.0, 0.0, 0.0]
    disp.LineWidth = 3.0   # thicker line reads as crisp black, not thin gray
    disp.SetScalarBarVisibility(view, False)

    # Decide 2-D vs 3-D from the data's actual bounding box (robust even if
    # you point this at an XDMF that isn't from writer.py)
    xmin, xmax, ymin, ymax, zmin, zmax = reader.GetDataInformation().GetBounds()
    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    is_3d = dz > 1e-9 * max(dx, dy, 1.0)

    cx, cy, cz = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0
    diag = (dx**2 + dy**2 + dz**2) ** 0.5 or 1.0
    d = diag * 2.0

    if is_3d:
        # Camera along the (1,1,1) diagonal from the data center = isometric
        view.CameraFocalPoint = [cx, cy, cz]
        view.CameraPosition   = [cx + d, cy + d, cz + d]
        view.CameraViewUp     = [0.0, 0.0, 1.0]
    else:
        # Flat mesh: look straight down the z-axis
        view.CameraFocalPoint = [cx, cy, cz]
        view.CameraPosition   = [cx, cy, cz + d]
        view.CameraViewUp     = [0.0, 1.0, 0.0]

    ResetCamera(view)   # keeps the direction/up just set, fits distance to data
    Render(view)

    SaveScreenshot(str(png_path), view, ImageResolution=[width, height])
    print(f"Saved {png_path} ({'3-D isometric' if is_3d else '2-D top-down'} view)")

    # Clean up so the next file starts from a blank view / fresh pipeline
    Hide(reader, view)
    Delete(disp)
    Delete(reader)


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: pvbatch render.py mesh.xdmf|folder/ [width] [height]")

    input_path = Path(sys.argv[1])
    width  = int(sys.argv[2]) if len(sys.argv) > 2 else 1600
    height = int(sys.argv[3]) if len(sys.argv) > 3 else 1200

    if input_path.is_dir():
        xdmf_files = sorted(input_path.glob("*.xdmf"))
        if not xdmf_files:
            sys.exit(f"No .xdmf files found in {input_path}")
    elif input_path.is_file():
        xdmf_files = [input_path]
    else:
        sys.exit(f"Not found: {input_path}")

    view = CreateView('RenderView')
    view.ViewSize = [width, height]
    view.Background = [1, 1, 1]
    view.OrientationAxesVisibility = 1
    configure_light_kit(view)

    for xdmf_path in xdmf_files:
        render_one(xdmf_path, view, width, height)


if __name__ == "__main__":
    main()