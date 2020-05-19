function[Treat]=ExtractTailH(frame)

Treat=frame;

[~, t0] = edge(Treat, 'sobel');
fudgeFactor = .7;
Treat = edge(Treat,'sobel', t0 * fudgeFactor);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
Treat = imdilate(Treat, [se90 se0]);
Treat = imfill(Treat, 'holes');
Treat = bwareaopen(Treat,500);

end


