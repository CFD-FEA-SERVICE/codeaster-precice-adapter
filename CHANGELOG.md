# Changelog

## 0.1.3
- **CHT fix (energy conservation at the interface).** `precice_aster.cht` now sends
  `Sink-Temperature-Solid = T_wall + (h_f / h_s) (T_wall - T_sink_fluid)` instead of the wall
  temperature. The OpenFOAM adapter reads the partner's sink temperature as the value at a
  distance delta inside the partner; sending the wall temperature added a spurious resistance
  delta/lambda on the fluid side only, so the fluid received less heat than the solid gave.
  On the heated-plate example the old version over-predicted the heat through the plate by
  about 2x and the interface came out 3-4 K too cold. With the fix the heat flux seen by the
  fluid equals the solid's for any `delta`.

## 0.1.2
- Code comments throughout `precice_aster`.
- New example `examples/elastic-tube-3d` (3D TETRA10, OpenFOAM `pimpleFoam`), with tightened
  fluid tolerances and absolute-or-relative convergence measures.
- New tool `tools/inp_to_mail.py`: Abaqus/CalculiX tetrahedral meshes to code_aster `.mail`.

## 0.1.1
- First public release: FSI (`DYNA_NON_LINE`) and CHT (`THER_LINEAIRE` / `THER_NON_LINE`).
