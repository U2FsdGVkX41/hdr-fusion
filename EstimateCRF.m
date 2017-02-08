function CRF = EstimateCRF(imageData, exposureTime)

zMax = 255;
zMin = 0;
zMid = (zMax + zMin) * 0.5;

weight = zeros(256,1);
for h = 1:256
    if (h-1) > zMid
        weight(h) = zMax - (h - 1);
    else
        weight(h) = (h - 1) - zMin;
    end
end

stack_hist = ComputeLDRStackHistogram(imageData);
stack_samples = GrossbergSampling('', '', stack_hist, 100);

intensityRange = 256;
CRF = zeros(intensityRange, 3);
for i=1:3
    CRF(:,i) = gsolve(stack_samples(:,:,i), log(exposureTime), 20, weight);
end

end