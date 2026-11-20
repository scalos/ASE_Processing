function tests = test_segToAseMembership
% Grid geometry of segToAseMembership: NIfTI affines map 0-based indices, MATLAB arrays are
% 1-based. With identical grids an index-convention error cancels, so these tests use grids that
% differ in voxel size, origin and axis direction.
    tests = functiontests(localfunctions);
end

function info = gridInfo(dims, vox, origin, flipY)
% niftiinfo-like struct: Transform.T (row-vector affine3d) maps 0-based [i j k 1] to world mm
    A = diag(vox);
    if flipY, A(2,2) = -A(2,2); end
    info.ImageSize = dims;
    info.Transform = affine3d([A zeros(3,1); origin(:)' 1]);
end

function [seg, LUT] = xLabelSeg(dims)
% label of each seg voxel = its 0-based x index + 1
    [I,~,~] = ndgrid(0:dims(1)-1, 0:dims(2)-1, 0:dims(3)-1);
    seg = I + 1;
    LUT = table((1:dims(1))', 'VariableNames', {'code'});
end

function lab = dominant(M, LUT, dims)
    [mx, r] = max(M, [], 1);
    lab = reshape(LUT.code(full(r)), dims);  lab(full(mx) == 0) = 0;
end

function testCentreSampleCoarseGrid(tc)
% ASE 2 mm, seg 1 mm, shared origin: ASE voxel i0 centre sits on seg voxel 2*i0 -> code 2*i0+1
    segDims = [20 6 6];  aseDims = [8 3 3];
    [seg, LUT] = xLabelSeg(segDims);
    si = gridInfo(segDims, [1 1 1], [0 0 0], false);
    ai = gridInfo(aseDims, [2 2 2], [0 0 0], false);
    M = segToAseMembership(seg, LUT, si, ai, 'maxSampling', [1 1 1]);
    lab = dominant(M, LUT, aseDims);
    i0 = (0:aseDims(1)-1)';
    verifyEqual(tc, lab(:,1,1), 2*i0 + 1);
end

function testOriginOffsetAndFlippedAxis(tc)
% ASE grid with a shifted origin and a reversed y axis; labels follow world x, so every ASE
% voxel's label is fixed by its world x: code = (origin_x + 2*i0) - segOrigin_x + 1
    segDims = [30 30 10];  aseDims = [10 10 3];
    [seg, LUT] = xLabelSeg(segDims);
    si = gridInfo(segDims, [1 1 1], [-5 -5 0], false);
    ai = gridInfo(aseDims, [2 2 3], [1 14 2], true);          % y runs downward from 14 mm
    M = segToAseMembership(seg, LUT, si, ai, 'maxSampling', [1 1 1]);
    lab = dominant(M, LUT, aseDims);
    i0 = (0:aseDims(1)-1)';
    expected = (1 + 2*i0) - (-5) + 1;
    verifyEqual(tc, lab(:,1,1), expected);
    verifyEqual(tc, lab(:,end,end), expected);                % independent of y and z
end

function testOccupancyHalfVoxel(tc)
% 2 mm ASE voxels over a 1 mm seg whose boundary splits every ASE voxel in half along x:
% each voxel is 50/50 between its two x-labels when sub-sampled
    segDims = [16 4 4];  aseDims = [8 2 2];
    [seg, LUT] = xLabelSeg(segDims);
    si = gridInfo(segDims, [1 1 1], [0 0 0], false);
    ai = gridInfo(aseDims, [2 2 2], [0.5 0.5 0.5], false);    % centres between seg voxels
    M = segToAseMembership(seg, LUT, si, ai, 'maxSampling', [2 2 2]);
    full1 = full(M(:, 1));
    verifyEqual(tc, nnz(full1), 2);
    verifyEqual(tc, full1(full1 > 0), [0.5; 0.5], 'AbsTol', 1e-12);
    verifyEqual(tc, full(sum(M, 1)), ones(1, prod(aseDims)), 'AbsTol', 1e-12);
end
