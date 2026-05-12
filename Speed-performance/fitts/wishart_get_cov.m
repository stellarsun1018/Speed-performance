function cov_mats = wishart_get_cov(phi,dist,dur)
wishart_sublayer = @(theta,dist,dur) theta(1) .* dist + theta(2) .* dur + theta(3);

relay_mat_flat = NaN(6,length(dist));
for i = 1:6
    relay_mat_flat(i,:) = wishart_sublayer(phi(i,:,:),dist,dur)';
end
relay_mat = reshape(relay_mat_flat,3,2,length(dist));

cov_mats = NaN(2,2,length(dist));
for i = 1:length(dist)
    cov_mats(:,:,i) = 0.1 .* squeeze(relay_mat(:,:,i))' * squeeze(relay_mat(:,:,i));
end

end