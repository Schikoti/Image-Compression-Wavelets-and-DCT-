clc;
clear all;
close all;
fid = fopen('lena.raw');
Y = fread(fid,[512,512],'uchar');
fclose(fid);
Y = Y';
rows=512;
columns=512;
for row=1:512
    for column=1:512
      X(row,column) = Y(row,column);
    end
end

  h1 = [0.026749, -0.016864,-0.078223,0.266864,0.602949,0.266864,-0.078223,-0.016864,0.026749];

  h2 = [0,-0.045636,-0.028772,0.295636,0.557543,0.295636,-0.028772,-0.045636,0];
  
  g1 = [-h2(9),h2(8),-h2(7),h2(6),-h2(5),h2(4),-h2(3),h2(2),-h2(1)];

  g2 = [-h1(9),h1(8),-h1(7),h1(6),-h1(5),h1(4),-h1(3),h1(2),-h1(1)];

  for n=1:9
    h1(n) = h1(n)*sqrt(2.0);
    h2(n) = h2(n)*sqrt(2.0);
    g1(n) = g1(n)*sqrt(2.0);
    g2(n) = g2(n)*sqrt(2.0);
  end
  %level1
  X_new=analysis_filter_rows(Y,512,512,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,512/2,512,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,512/2,512,256,1,h1,g1);
  %level2 top left
  X_new=analysis_filter_rows(X_new,256,256,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,128,1,h1,g1);
  %level2 top right
  X_new=analysis_filter_rows(X_new,256,256,1,256,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,1,256,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,128,256,h1,g1);
  %level2 lower left
  X_new=analysis_filter_rows(X_new,256,256,256,1,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,256,1,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,128+256,1,h1,g1);
  %level2 lower right
  X_new=analysis_filter_rows(X_new,256,256,256,256,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,256,256,h1,g1);
  X_new=analysis_filter_columns(X_new,128,256,128+256,256,h1,g1);
  %level3 divide lower left-hand subband into 4 subbands creating 19 subbands
  X_new=analysis_filter_rows(X_new,128,128,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,64,128,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,64,128,64,1,h1,g1);
  % level4 Further divide lower left-hand subband into 4 subbands creating 22 subbands
  X_new=analysis_filter_rows(X_new,64,64,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,32,64,1,1,h1,g1);
  X_new=analysis_filter_columns(X_new,32,64,32,1,h1,g1);
  
  X_new_3=zeroed_22(X_new,3);
  X_new_10=zeroed_22(X_new,10);
  X_new_15=zeroed_22(X_new,15);
  X_new_3=quantization_22band(X_new_3);
  X_new_10=quantization_22band(X_new_10);
  X_new_15=quantization_22band(X_new_15);
  X_new_3=zeroed_22(X_new_3,3);
  X_new_10=zeroed_22(X_new_10,10);
  X_new_15=zeroed_22(X_new_15,15);
%   
%   % synthesis level4
%   

Y_rec=synthesis_filter_columns(X_new_3,32,32,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/16,columns/16,rows/16,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/16,columns/8,1,1,h2,g2);
  %level3
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,rows/8,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/8,columns/4,1,1,h2,g2);
  %level2 upper left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,1,h2,g2);
  %level2 upper right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,columns/2,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,columns/2,h2,g2);
  %level2 lower left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,1,h2,g2);
  %level2 lower right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,columns/2,h2,g2);  
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,columns/2,h2,g2);
  %level1
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/2,columns,1,1,h2,g2);
  for row=1:rows
      for column=1:columns
          if Y_rec(row,column)<0
              Y_rec(row,column)=0;
          end
          if Y_rec(row,column)>255
              Y_rec(row,column)=255;
          end
      end
  end
  psnr3=calc_psnr(Y_rec,X);
  subplot(2,4,1)
  imshow(X,[])
  title('Original Lina Image')
  subplot(2,4,2)
  imshow(Y_rec,[])
  title('Reconstructed Image1(Three-highest freq-Zero)')
  
  Y_rec=synthesis_filter_columns(X_new_10,32,32,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/16,columns/16,rows/16,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/16,columns/8,1,1,h2,g2);
  %level3
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,rows/8,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/8,columns/4,1,1,h2,g2);
  %level2 upper left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,1,h2,g2);
  %level2 upper right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,columns/2,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,columns/2,h2,g2);
  %level2 lower left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,1,h2,g2);
  %level2 lower right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,columns/2,h2,g2);  
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,columns/2,h2,g2);
  %level1
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/2,columns,1,1,h2,g2);
  for row=1:rows
      for column=1:columns
          if Y_rec(row,column)<0
              Y_rec(row,column)=0;
          end
          if Y_rec(row,column)>255
              Y_rec(row,column)=255;
          end
      end
  end
  psnr10=calc_psnr(Y_rec,X);
  subplot(2,4,3)
  imshow(Y_rec,[])
  title('Reconstructed Image2(Ten-highest freq-Zero)')
  
  Y_rec=synthesis_filter_columns(X_new_15,32,32,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/16,columns/16,rows/16,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/16,columns/8,1,1,h2,g2);
  %level3
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/8,columns/8,rows/8,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/8,columns/4,1,1,h2,g2);
  %level2 upper left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,1,h2,g2);
  %level2 upper right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,1,columns/2,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,1,columns/2,h2,g2);
  %level2 lower left
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,1,h2,g2);
  %level2 lower right
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2,columns/2,h2,g2);  
  Y_rec=synthesis_filter_columns(Y_rec,rows/4,columns/4,rows/2 + rows/4,columns/2,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/4,columns/2,rows/2,columns/2,h2,g2);
  %level1
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,1,1,h2,g2);
  Y_rec=synthesis_filter_columns(Y_rec,rows/2,columns/2,rows/2,1,h2,g2);
  Y_rec=synthesis_filter_rows(Y_rec,rows/2,columns,1,1,h2,g2);
  for row=1:rows
      for column=1:columns
          if Y_rec(row,column)<0
              Y_rec(row,column)=0;
          end
          if Y_rec(row,column)>255
              Y_rec(row,column)=255;
          end
      end
  end
  psnr15=calc_psnr(Y_rec,X);
  subplot(2,4,4)
  imshow(Y_rec,[])
  title('Reconstructed Image3(fifteen-highest freq-Zero)')
  
%   Magnitude spectrum
X_fft=fft_basic(fft_basic(X).').';
X_fft_3=fft_basic(fft_basic(X_new_3).').';
X_fft_10=fft_basic(fft_basic(X_new_10).').';
X_fft_15=fft_basic(fft_basic(X_new_15).').';
subplot(2,4,5)
imshow(log(fftshift(abs(X_fft))),[])
title('2D-DFT Magnitude Spectrum of Original Lina Image')
subplot(2,4,6)
imshow(log(fftshift(abs(X_fft_3))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image1')
subplot(2,4,7)
imshow(log(fftshift(abs(X_fft_10))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image2')
subplot(2,4,8)
imshow(log(fftshift(abs(X_fft_15))),[])
title('2D-DFT Magnitude Spectrum of Reconstructed image3')

function X_new= zeroed_22(X_new, num)
if num==3
    for i=257:384
        for j=385:512
            X_new(i,j)=0;
        end
    end
    for i=385:512
        for j=257:512
            X_new(i,j)=0;
        end
    end
end
if num==10
    for i=1:512
        for j=385:512
            X_new(i,j)=0;
        end
    end
    for i=129:512
        for j=257:384
            X_new(i,j)=0;
        end
    end
    for i=257:512
        for j=129:256
            X_new(i,j)=0;
        end
    end
    for i=385:512
        for j=1:128
            X_new(i,j)=0;
        end
    end
end

if num==15
    for i=1:512
        for j=129:512
            X_new(i,j)=0;
        end
    end
    for i=129:512
        for j=1:128
            X_new(i,j)=0;
        end
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

function psnr=calc_psnr(image,lina)
mse=sum(sum((image-lina).^2))/(512*512);
psnr=10*log10(255*255/mse);
end

function quant_22 = quantization_22band(ac)

minimum_ac = min(min(ac));

step_ac = max((max(ac))-min(min(ac)))/1024;

for i=1:512
    for j=1:512
        temp_ac = floor((ac(i,j)-minimum_ac)/step_ac);
        quant_22(i,j)=((temp_ac*step_ac+minimum_ac)+(temp_ac*step_ac+minimum_ac))/2;
    end
end
end

function or_img= synthesis_filter_rows(or_img,rows,columns,row_offset,column_offset,h,g)
N = rows*2;
delay=5;%might change
temp=[];
filter_taps=9;
for column=column_offset:(columns+column_offset-1)
    for n=1:delay-1
        x1(n) = or_img((delay-1-n+row_offset),column);
    end
    for n=1:rows
        x1(n+delay-1) = or_img((n+row_offset-1),column);
    end
    
    for n=rows+1:((rows+delay-1))
        x1(n+delay-1) = or_img((2*rows-n-1+row_offset),column);
    end

    for n=1:rows+2*(delay-1)
          temp(2*n) = 0;
          temp(2*n-1) = x1(n);
    end
 
    for m=1:2:N  
        sum = 0;
        for n=(m-filter_taps):2:m-1
            sum= sum+temp(n+filter_taps)*h(m-n);
        end
        z1(m) = sum;
        sum = 0;
        for n=((m+2)-filter_taps):2:m-1
            sum= sum+temp(n+filter_taps)*h((m+1)-n);
        end
        z1(m+1) = sum;
    end
 
    for n=1:delay-1
        x2(n) = or_img((delay-1-n+row_offset+rows),column);
    end

    for n=1:rows
        x2(n+delay-1) = or_img((n+row_offset+rows-1),column);
    end

    for n=rows+1:(rows+delay-1)
        x2(n+delay-1) = or_img((2*rows-n-1+row_offset+rows),column);
    end
 
    for n=1:((rows+2*(delay-1)))
        temp(2*n-1) = 0;
        temp(2*n) = x2(n);    
    end
 
    for m=1:2:N
        sum=0;
        for n=((m+1)-filter_taps):2:m-1
            sum=sum+temp(n+filter_taps)*g(m-n);
        end
        z2(m)=sum;

        sum = 0;
        for n=((m+1)-filter_taps):2:(m)
            sum=sum+temp(n+filter_taps)*g((m+1)-n);
        end

        z2(m+1) = sum;

    end
 
    for row=1:N
          or_img((row+row_offset-1),column) = z1(row)+ z2(row);
    end
end
%     for n=1:ceil(delay/2)
%         x1(n) = or_img((ceil(delay/2)-n+row_offset+2),column);
%     end
%     for n=1:rows
%         x1(n+ceil(delay/2)) = or_img((n+row_offset-1),column);
%     end
%     
%     for n=rows+1:((rows+ceil(delay/2)))
%         x1(n+ceil(delay/2)) = or_img((2*rows-n-1+row_offset),column);
%     end
% 
%     for n=1:rows+delay+1
%           temp(2*n) = 0;
%           temp(2*n-1) = x1(n);
%     end
%  
%     for m=1:2:N  
%         sum = 0;
%         for n=(m-filter_taps):2:m-1
%             sum= sum+temp(n+filter_taps)*h(m-n);
%         end
%         z1(m) = sum;
%         sum = 0;
%         for n=((m+2)-filter_taps):2:m-1
%             sum= sum+temp(n+filter_taps)*h((m+1)-n);
%         end
%         z1(m+1) = sum;
%     end
%  
%     for n=1:ceil(delay/2)
%         x2(n) = or_img((ceil(delay/2)-n+2+row_offset+rows),column);
%     end
% 
%     for n=1:rows
%         x2(n+ceil(delay/2)) = or_img((n+row_offset+rows-1),column);
%     end
% 
%     for n=rows+1:(rows+ceil(delay/2))
%         x2(n+ceil(delay/2)) = or_img((2*rows-n-1+row_offset+rows),column);
%     end
%  
%     for n=1:((rows+delay+1))
%         temp(2*n-1) = 0;
%         temp(2*n) = x2(n);    
%     end
%  
%     for m=1:2:N
%         sum=0;
%         for n=((m+1)-filter_taps):2:m-1
%             sum=sum+temp(n+filter_taps)*g(m-n);
%         end
%         z2(m)=sum;
% 
%         sum = 0;
%         for n=((m+1)-filter_taps):2:(m)
%             sum=sum+temp(n+filter_taps)*g((m+1)-n);
%         end
% 
%         z2(m+1) = sum;
% 
%     end
%  
%     for row=1:N
%           or_img((row+row_offset-1),column) = z1(row)+ z2(row);
%     end
% end
end

 function or_img=synthesis_filter_columns(or_img,rows,columns,row_offset,column_offset,h,g)

N = columns*2;
delay=5;
filter_taps=9;
for row=row_offset:(rows+row_offset)-1
    for n=1:delay-1
        x1(n) = or_img(row, ((delay-1)-n+column_offset));
    end
    for n=1:columns
        x1(n+(delay-1)) = or_img(row,(n+column_offset-1));
    end
    for n=columns+1:(columns+(delay-1))
        x1(n+(delay-1)) = or_img(row, ((2*columns)-n-1+column_offset));
    end
    for n=1:(columns+2*(delay-1))
        temp(2*n) = 0;
        temp(2*n-1) = x1(n);
    end   
    for m=1:2:N
        sum=0;
        for n=(m-filter_taps):2:m-1
            sum=sum+temp(n+filter_taps)*h(m-n);
        end
        z1(m)=sum;
        sum=0;
        for n=(m+2-filter_taps):2:m-1
            sum=sum+temp(n+filter_taps)*h(m+1-n);
        end
        z1(m+1)=sum;
    end
    for n=1:delay-1
        x2(n) = or_img(row, ((delay-1)-n+column_offset+columns));
    end
    for n=1:columns
        x2(n+(delay-1)) = or_img(row,(n+column_offset-1+columns));
    end
    for n=columns+1:(columns+(delay-1))
        x2(n+(delay-1)) = or_img(row, ((2*columns)-n-1+column_offset+columns));
    end
    for n=1:(columns+2*(delay-1))
        temp(2*n-1) = 0;
        temp(2*n) = x2(n);
    end 
    for m=1:2:N
        sum=0;
        for n=(m+1-filter_taps):2:m-1
            sum=sum+temp(n+filter_taps)*g(m-n);
        end
        z2(m)=sum;
        sum=0;
        for n=(m+1-filter_taps):2:m-1
            sum=sum+temp(n+filter_taps)*g(m+1-n);
        end
        z2(m+1)=sum;
    end
    for column=1:N
       or_img(row, (column+column_offset-1)) = z1(column) + z2(column);
     end
     end
% for row=row_offset:(rows+row_offset-1)  
%     for n=1:ceil(delay/2)
%         x1(n) = or_img(row,(ceil(delay/2)+2-n+column_offset));
%     end
%     for n=1:columns
%         x1(n+ceil(delay/2)) = or_img(row,(n+column_offset-1));
%     end
%     for n=columns+1:(columns+ceil(delay/2))
%         x1(n+ceil(delay/2)) = or_img(row, (2*columns-n-1+column_offset));
%     end
%     for n=1:(columns+delay+1)
%         temp(2*n) = 0;
%         temp(2*n-1) = x1(n);
%     end
%     for m=1:2:N 
%         sum = 0;
%         for n=(m-filter_taps):2:m-1
%             sum= sum+temp(n+filter_taps)*h(m-n);
%         end
%         z1(m) = sum;
%         sum = 0;
%         for n=((m+2)-filter_taps):2:m-1
%             sum= sum+temp(n+filter_taps)*h((m+1)-n);
%         end
%             z1(m+1) = sum;
%     end
% 
%     for n=1:ceil(delay/2)
%         x2(n) = or_img(row,(ceil(delay/2)-n+2+column_offset+columns));
%     for n=1:columns
%         x2(n+ceil(delay/2)) = or_img(row,(n+column_offset+columns-1));
%     end
%     for n=columns+1:(columns+ceil(delay/2))
%       x2(n+ceil(delay/2)) = or_img(row,((2*columns)-n-1+column_offset+columns));
%     end
%     for n=1:(columns+delay+1)
%       temp(2*n-1) = 0;
%       temp(2*n) = x2(n);
%     end
%     for m=1:2:N
%       sum = 0;
%       for n=((m+1)-filter_taps):2:m-1
%         sum= sum+temp(n+filter_taps)*g(m-n);
%       end
%       z2(m) = sum;
%       sum = 0;
%       for n=((m+1)-filter_taps):2:(m)
%         sum=sum+temp(n+filter_taps)*g((m+1)-n);
%       end
%       z2(m+1) = sum;
%     end
% 
%     for column=1:N
%       or_img(row, (column+column_offset-1)) = z1(column) + z2(column);
%     end
%     end
% end

end

function or_img=analysis_filter_rows(or_img,rows,columns,row_offset,column_offset,h,g)
N = rows;
delay=5;
filter_taps=9;
for column=column_offset:columns+column_offset-1
    for n=1:delay-1
      x1(n) = or_img(((delay-1)-n+row_offset),column);
    end
    
    for n=1:N
        x1(n+(delay-1)) = or_img((n+row_offset-1),column);
    end
    
    for n=N+1:(N+(delay-1))
      x1(n+(delay-1)) = or_img(((2*N)-n-1+row_offset),column);
    end
    
    for n=1:(N+2*(delay-1))
      x2(n) = x1(n);
    end

    for m=1:2:N
      sum = 0;
      for n=(m-filter_taps):m-1
        sum= sum+(x1(n+filter_taps))*h(m-n);
      end
      z1(m) = sum;
    end

    for m=2:2:N
      sum = 0;
      for n=(m-filter_taps):m-1
        sum= sum+(x2(n+filter_taps))*g(m-n);
      end
      z2(m) = sum;
     end

 %decimate and write back low-frequency result%
    for row=1:(rows/2)
      or_img((row+row_offset-1),column) = z1(row*2-1);
    end

%decimate and write back high-frequency result%
    for row=1:(rows/2)
      or_img((row+row_offset+rows/2-1),column) = z2(row*2);
    end
end
end


function or_img=analysis_filter_columns(or_img,rows,columns,row_offset,column_offset,h,g)

N = columns;
delay=5;
filter_taps = 9;
for row=row_offset:(rows+row_offset)-1
    for n=1:delay-1
        x1(n) = or_img(row, ((delay-1)-n+column_offset));
    end
    for n=1:N
        x1(n+(delay-1)) = or_img(row,(n+column_offset-1));
    end
    for n=N+1:(N+(delay-1))
        x1(n+(delay-1)) = or_img(row, ((2*N)-n-1+column_offset));
    end
    for n=1:(N+2*(delay-1))
        x2(n) = x1(n);
    end
    for m=1:2:N
        sum = 0;
        for n=(m-filter_taps):m-1
            sum= sum+x1(n+filter_taps)*h(m-n);
        end
        z1(m) = sum;
    end

    for m=2:2:N  
      sum = 0;
      for n=(m-filter_taps):m-1
        sum= sum+x2(n+filter_taps)*g(m-n);
      end
      z2(m) = sum;
    end
    

    %decimate and write back low-frequency result%
    
    for column=1:(columns/2)
      or_img(row, (column+column_offset-1)) = z1(column*2-1);
    end


    %decimate and write back high-frequency result%
    for column=1:(columns/2)
      or_img(row,(column+column_offset+columns/2-1)) = z2(column*2);
    end
end



end
