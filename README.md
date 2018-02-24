## Dependency

- cmake
- eigen3
- glog
- gflags
- tbb
- boost
- opencv 3.0
- mpfr
- gomp
- suitesparse
- ceres-solver
- qt4
- glew

## Usage 

1. Download data from [project website](http://www.yf.io/p/tiny). Assume the data is downloaded to `<tiny_dir>/images/`.
2. Detect correspondences with KLT

```
./tiny_detect_tracks -in_image_dir <tiny_dir>/images/stone6_still/ -in_image_pattern ".*" \
    -num_images 100 -logtostderr=1 -max_corners 2000 -corner_block_size 5  -corner_min_distance 3 \
    -out_model <tiny_dir>/models/stone6_klt -noshow -lk_error_threshold 6
```

2. Bundle adjustment to get camera poses and structures
```
./tiny_bundle_adjust -in_model <tiny_dir>/models/stone6_klt -out_model <tiny_dir>/models/stone6_ba1 \
    -settings settings/ba_settings.txt -set_p_random -set_c_zero --v=1
```

3. Create multi-view stereo cost volume

```
./tiny_sweep_planes -in_model <tiny_dir>/models/stone6_ba1 -scale 2 -patch_radius 1 -num_samples 64 \
    -out_cv <tiny_dir>/volumes/stone6_s2_p3_n_64.cv -logtostderr
```

4. Use CRF to get smooth depth estimation
```
./tiny_solve_depth -in_cv <tiny_dir>/volumes/stone6_s2_p3_n_64.cv -out_dm stone6_drm.dm -logtostderr \
    -out_image <tiny_dir>/output/stone6_crf.png -use_densecrf -settings settings/dense_crf_settings.txt \
    -smooth_settings settings/smooth_settings.txt
```
