% 自检5/6 失配诊断(面元位置修复后)
addpath('scene','signal');
c = 299792458;
cfg = config(); rng(cfg.sim.seed);
lambda = c/cfg.waveform.fc;
wf = Waveform_derive(cfg);
sc.radar(1).p_true = [0 0 0];  sc.radar(1).R_true = eye(3);
sc.radar(1).boresight = [1;0;0];
sc.radar(1).ntx = 1;  sc.radar(1).nrx = 1;
sc.radar(1).d_tx = lambda/2;  sc.radar(1).d_rx = lambda/2;
sc.radar(1).tx_pos = [0;0;0];  sc.radar(1).rx_pos = [0;0;lambda/2];
sc.radar(1).fov = cfg.radar(1).fov;  sc.radar(1).g_p = cfg.radar(1).g_p;
Pt = 10^(cfg.radar(1).Pt_dbm/10)*1e-3;
G0 = 10^(cfg.radar(1).G0_dbi/10);
Pn = 10^((wf.noise_dbm-30)/10);
N  = wf.N_fast;
Pa_h = zeros(1,2);
for cube = [1 2]
    cfg5 = cfg;
    cfg5.target.targets(1).pos  = [4.5 + (cube-1)*0.5, 0, 0];
    cfg5.target.targets(1).size = [cube cube cube];
    tgt = Scene_targets(cfg5);
    e  = Echo_generate(cfg5, sc, tgt, wf, true);
    s1 = e(1).data(:,1,1);
    X  = fft(s1);
    Pa = sum(abs(X).^2);
    Pa_h(cube) = Pa;
    St = N*sum(abs(s1).^2);
    Pm = tgt(1).facet.pos.';   Nm = tgt(1).facet.nrm.';   Af = tgt(1).facet.area;
    r   = vecnorm(Pm,2,1).';
    uTx = -Pm./r.';
    nTx = sum(Nm.*uTx, 1).';
    cosb = max(nTx,0);
    sig = 4*pi/lambda^2 * Af.^2 .* cosb.^2;
    Prx = Pt*G0^2*lambda^2*sig./((4*pi)^3 * r.^2 .* r.^2);
    a = sqrt(max(Prx,0)/Pn);   a(cosb<=0) = 0;
    tau = 2*r/c;
    phi = 2*pi*(wf.fc*tau - wf.S/2.*tau.^2);
    S1 = sum(a.*exp(1i*phi));
    asum = sum(a.^2);
    fprintf('== cube %g m: M=%d\n', cube, tgt(1).M);
    fprintf('  Pa(fft)=%.1f dB | St(N*sum|s|^2)=%.1f dB\n', 10*log10(Pa), 10*log10(St));
    fprintf('  N^2*sum(a^2)=%.1f dB | N^2*|S1|^2=%.1f dB | asum=%.4g\n', ...
            10*log10(N^2*asum), 10*log10(N^2*abs(S1)^2), asum);
    [pk, ki] = sort(abs(X).^2, 'descend');
    fprintf('  top6 bins(0-based:power_dB):');
    for q = 1:6
        fprintf(' [%d:%.1f]', ki(q)-1, 10*log10(pk(q)));
    end
    fprintf('\n');
end
fprintf('Pb/Pa = %.2f dB\n', 10*log10(Pa_h(2)/Pa_h(1)));
