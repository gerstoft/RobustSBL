function level = maximum_sidelobe_level(b,phi_vec,phi_src,phi_tol)
%level = maximum_sidelobe_level(b,phi_vec,phi_src,phi_toler)
%
%evaluation of maximum sidelobe level in beampattern b
%
%INPUT
%b        is a vector with samples of the beampattern
%phi_vec  is a vector with bearings, same size as b
%phi_src  is the true direction of arrival, a scalar value
%phi_tol  is the tolerance in direction finding
%
%OUTPUT
%level    is a scalar, the estimated maximum sidelobe level

[pks,locs]     = findpeaks(b,phi_vec);

[level,index] = max(pks); % find highest peak in beampattern b

% check if the highest peak is near the true source DOA

if abs(locs(index) - phi_src) < phi_tol,
    pks(index) = []; % delete the correct main lobe
    locs(index)= []; % delete the valid DOA of the main lobe
    level = max(pks); % find second highest peak in beampattern b
else
    % the highest peak is not near the true source DOA 
    % in this case, we return the maximum peak of the beampattern b
end

if isempty(level)
    level = 0;
end
