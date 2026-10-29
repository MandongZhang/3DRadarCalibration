function P = Det_to_points(D)
% 检测/真值记录 (r, az, el) -> 各站站系直角坐标点(点云累积用)
%
% 输入: D (Nf, n_st, K>=3),前三切片 = [r (m) | az (deg) | el (deg)],站系球坐标
% 输出: P (Nf, n_st, 3) 站系直角坐标,与 Truth_project 的投影互逆:
%   x = r*cos(el)*cos(az),  y = r*cos(el)*sin(az),  z = r*sin(el)

    Nf  = size(D, 1);
    n_st = size(D, 2);
    r   = D(:, :, 1);
    az  = D(:, :, 2);
    el  = D(:, :, 3);
    P(:, :, 1) = r .* cosd(el) .* cosd(az);
    P(:, :, 2) = r .* cosd(el) .* sind(az);
    P(:, :, 3) = r .* sind(el);
end
