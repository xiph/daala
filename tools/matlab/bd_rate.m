#!/usr/bin/octave -qf

warning("off","Octave:nested-functions-coerced");

args=argv();

if size(args,1)!=2
  printf("usage: ./bd_rate.sh <RD-1.out> <RD-2.out>\n");
  return
end

TYPE=getenv("TYPE");
if strcmp(TYPE,"")
  TYPE="piecewise-linear";
end

switch (TYPE)
  case "piecewise-linear"
    t=1;
  case "cubic-polyfit"
    t=2;
%  case "two-point-monotone-cubic-spline"
%    t=3;
%  case "three-point-monotone-cubic-spline"
%    t=4;
%  case "cubic-spline-interp"
%    t=5;
%  case "shape-preserving-cubic-hermite-interp"
%    t=6;
  otherwise
    printf("Invalid type: %s\n",TYPE);
    return
endswitch

rd1=load("-ascii",args{1});
rd2=load("-ascii",args{2});

rd1=flipud(sortrows(rd1,1));
rd2=flipud(sortrows(rd2,1));

rate1=rd1(:,3)*8./rd1(:,2);
rate2=rd2(:,3)*8./rd2(:,2);

pin = program_invocation_name;
chdir(pin(1:(length(pin)-length(program_name))));

[psnr_rate,psnr_dsnr]=bjontegaard([rate1,rd1(:,4)],[rate2,rd2(:,4)],t);
[psnrhvs_rate,psnrhvs_dsnr]=bjontegaard([rate1,rd1(:,5)],[rate2,rd2(:,5)],t);
[ssim_rate,ssim_dsnr]=bjontegaard([rate1,rd1(:,6)],[rate2,rd2(:,6)],t);
[fastssim_rate,fastssim_dsnr]=bjontegaard([rate1,rd1(:,7)],[rate2,rd2(:,7)],t);

printf("           RATE (%%)  DSNR (dB)\n");
printf("    PSNR %0.5f  %0.5f\n",psnr_rate,psnr_dsnr);
printf(" PSNRHVS %0.5f  %0.5f\n",psnrhvs_rate,psnrhvs_dsnr);
printf("    SSIM %0.5f  %0.5f\n",ssim_rate,ssim_dsnr);
printf("FASTSSIM %0.5f  %0.5f\n",fastssim_rate,fastssim_dsnr);
