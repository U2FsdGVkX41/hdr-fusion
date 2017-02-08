function pixel_value = iCRF(CRF, radiance, T)

rad_column = radiance(:);

l = log(rad_column*T);
diff = zeros(length(rad_column), 256);

for n=1:256
    diff(:,n) = abs(l-CRF(n));
end 

[~,pixel_value] = min(diff, [], 2);

pixel_value = pixel_value - 1;

pixel_value = reshape(pixel_value, size(radiance));
end

