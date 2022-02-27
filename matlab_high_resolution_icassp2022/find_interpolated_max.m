function doa = find_interpolated_max(phi_vec,b,index)
doa = phi_vec(index); % doa for maximum sample
if (index > 1) && (index < length(phi_vec)),
    aux = index + [-1 0 1];
    [~, doa] = parabel(phi_vec(aux),b(aux)); % doa for interpolated maximum
end
