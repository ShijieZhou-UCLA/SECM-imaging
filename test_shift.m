figure();

for i = 1:16
    subplot(4,4,i)
    
    bplines   = ScanLines(lines);
    bplines.currents = circshift(bplines.currents, -i*10, 1);

    p0        = ProbeParams(0);
    p0.angles = ProbeParam(angles);
    % p0.centerbias.set_value(lines.params.centerbias.value)
    bplines.params = p0;
    bpimage = bplines.back_project();
    bpimage.draw_image();
end
