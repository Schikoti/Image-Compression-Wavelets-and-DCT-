% IMAGE COMPRESSION using Discreate Cosine Transform
% Given Input image has been extracted as a 2D array for processing it further. 

clc;
clear all;
close all;
fid = fopen('lena.raw');
lina_array = fread(fid,[512,512],'uchar');
fclose(fid);
lina_array = lina_array';


% Initializing the array's.
%calculating the DCT blockwise using the function dct_test which is written from scratch.
% Making the quarter, half and threequarters of DCT matrix to zero using
% the function makezeros.
[n,l]=size(lina_array);
dct_quarter=zeros(512);
dct_half=zeros(512);
dct_threequarter=zeros(512);
dct_array=zeros(512);
for i=1:8:n
    for j=1:8:l
        lina_block=lina_array(i:i+7,j:j+7);% accessing 8x8 blocks of Lina array.
        dct_array(i:i+7,j:j+7)=dct_basic(lina_block,0);% Obtaining the DCT matrix of lina image using dct developed from scratch.
        dct_quarter(i:i+7,j:j+7)=makezeros(dct_array(i:i+7,j:j+7),0);% Setting the quarter of high frequencies to zero using makezeros function.
        dct_half(i:i+7,j:j+7)=makezeros(dct_array(i:i+7,j:j+7),1);%Setting the half of high frequencies to zero.
        dct_threequarter(i:i+7,j:j+7)=makezeros(dct_array(i:i+7,j:j+7),2);%Setting the threequarters of high frequencies to zero.
    end
end

% extracting DC(Lowest frequency) values from the DCT 8x8 blocks.
dct_orig_low=zeros(64);
a=1;b=1;
for i=1:8:n
    b=1;
    for j=1:8:l
        dct_orig_low(a,b)=dct_array(i,j);% DC Values have been extracted from each 8x8 blocks
        b=b+1;
    end
    a=a+1;
end

% Quantizing the quarter, half and three quarter DCT arrays to 10 bit using
% the dct_quantization function.
dct_quarter_quantised=dct_quantization(dct_orig_low,dct_quarter);
dct_half_quantised=dct_quantization(dct_orig_low,dct_half);
dct_threequarter_quantised=dct_quantization(dct_orig_low,dct_threequarter);

% Quantization has been performed on the dct-quaterzeros, dct-halfzeros and dct-threequarterzeros matrices.
% Due to quantization the zeroed elements will also be replaced with quantized values.
%Hence quarter, half and threequarter portion of the quantized matrices have been set to zero again.
for i=1:8:n
    for j=1:8:l
        dct_quarter_QMZ(i:i+7,j:j+7)=makezeros(dct_quarter_quantised(i:i+7,j:j+7),0);
        dct_half_QMZ(i:i+7,j:j+7)=makezeros(dct_half_quantised(i:i+7,j:j+7),1);
        dct_threequarter_QMZ(i:i+7,j:j+7)=makezeros(dct_threequarter_quantised(i:i+7,j:j+7),2);
    end
end

% Reconstruction of the image from the compressed forms(Quarter, half and three-Quarter retained portions).
% IDCT of quarter, half and threequarter(quantized and made zeros(QMZ)) 
% Inverse DCT has also been calculated from scratch as a part dct_basic function. 
for i=1:8:n
    for j=1:8:l
        lina1QMZ(i:i+7,j:j+7)=dct_basic(dct_quarter_QMZ(i:i+7,j:j+7),1);
        lina2QMZ(i:i+7,j:j+7)=dct_basic(dct_half_QMZ(i:i+7,j:j+7),1);
        lina3QMZ(i:i+7,j:j+7)=dct_basic(dct_threequarter_QMZ(i:i+7,j:j+7),1);
    end
end

% PSNR Signal to noise ratio has been calculated between the original image
% and the three reconstructed images.
psnr1_QMZ=calc_psnr(lina1QMZ,lina_array);% calc_psnr funtion has been defined with the basic SNR formula and returns the PSNR.
psnr2_QMZ=calc_psnr(lina2QMZ,lina_array);
psnr3_QMZ=calc_psnr(lina3QMZ,lina_array);

% Plots of Original and reconstructed images and their 2D-DCT magnitude spectrums.
lina_fft=fft_basic(fft_basic(lina_array).').';
lina1QMZ_fft=fft_basic(fft_basic(lina1QMZ).').';
lina2QMZ_fft=fft_basic(fft_basic(lina2QMZ).').';
lina3QMZ_fft=fft_basic(fft_basic(lina3QMZ).').';
figure(1)
subplot(2,4,1)
imshow(lina_array,[])
title('Original Lina Image')
subplot(2,4,2)
imshow(lina1QMZ,[])
title('Reconstructed Lina image1-(One-Quarter-zeros)')
subplot(2,4,3) 
imshow(lina2QMZ,[])
title('Reconstructed Lina image2-(One-Half-zeros)')
subplot(2,4,4)
imshow(lina3QMZ,[])
title('Reconstructed Lina image3-(ThreeQuarter-zeros)')
subplot(2,4,5)
imshow(log(fftshift(abs(lina_fft))),[])
title('2D-DFT Magnitude Spectrum of original image')
subplot(2,4,6)
imshow(log(fftshift(abs(lina1QMZ_fft))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image1')
subplot(2,4,7)
imshow(log(fftshift(abs(lina2QMZ_fft))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image2')
subplot(2,4,8)
imshow(log(fftshift(abs(lina3QMZ_fft))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image3')


%labels for plots

% inputs and outputs in each function
% explanation to each step
% explanation of what is being done in each funtion 

function psnr=calc_psnr(image,lina)
mse=sum(sum((image-lina).^2))/(512*512);
psnr=10*log10(255*255/mse);
end

function dct_2d=dct_basic(lina_array,num)
N=8;
y=zeros(N,N);
n=0:7;
for k=0:N-1
    y(k+1,:)=cos((2*(n)+1)*pi*k/(2*N));
end
y(1,:)=y(1,:)/sqrt(N);
y(2:N,:)=y(2:N,:)*sqrt(2/N);
if num == 0
    dct_2d=y*lina_array*y';
else 
    dct_2d=y'*lina_array*y;
end

end

function dct_qt2 = dct_quantization(dc,ac)
minimum_dc = min(min(dc));
minimum_ac = min(min(ac));
step_dc = max((max(dc))-min(min(dc)))/1024;
step_ac = max((max(ac))-min(min(ac)))/1024;
for i=1:64
    for j=1:64
        temp_dc = floor((dc(i,j)-minimum_dc)/step_dc);
        dct_qt1(i,j) = ((temp_dc*step_dc+minimum_dc)+(temp_dc*step_dc+minimum_dc))/2;
    end
end
for i=1:512
    for j=1:512
        temp_ac = floor((ac(i,j)-minimum_ac)/step_ac);
        dct_qt2(i,j)=((temp_ac*step_ac+minimum_ac)+(temp_ac*step_ac+minimum_ac))/2;
    end
end
a=1;
for i=1:8:512
    b=1;
    for j=1:8:512
        dct_qt2(i,j)=dct_qt1(a,b);
        b=b+1;
    end
    a=a+1;
end
end

function matt=makezeros(matt,num)
ind=reshape(1:numel(matt), size(matt));
ind = fliplr( spdiags( fliplr(ind) ) ); 
ind(:,1:2:end) = flipud( ind(:,1:2:end) );
ind(ind==0) = [];                           
[~,n]=size(ind);
if num==0
    f=48;
end
if num==1
    f=32;
end
if num==2
    f=16;
end
for z=1:n
    if z>f
        matt(ind(z))=0;
    end
end
end

function [op] = fft_basic(ip)
N=length(ip);
n = 0:1:N-1; % row vector for n
k = 0:1:N-1; % row vecor for k
WN = exp(-1j*2*pi/N); % Twiddle factor (w)
nk = n'*k; % creates a N by N matrix of nk values
W = WN .^ nk; % DFT matrix
op = (W*ip );
end