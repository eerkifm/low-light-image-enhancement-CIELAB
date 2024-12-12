%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                                                                       
%   Proposed Method - Image enhancement method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

fprintf ('Proposed Method - Linear Transformation Image Enhancement Method Based on CIELAB Color Space \n');
fprintf ('=================================================================\n\n');

%% ORIGINAL IMAGE TO XYZ WITH GAMMA CORRECTION
ori = imread('lowlight.bmp');      % Upload image
[N, M, P] = size (ori);     % Image size

ori_R = ori(:,:,1);       % Red plane
ori_G = ori(:,:,2);       % Green plane
ori_B = ori(:,:,3);       % Blue plane

% change original RGB to double
R_double = double(ori_R);
G_double = double(ori_G);
B_double = double(ori_B);

% ===========================================================================
% Reverse gamma correction to the original RGB
% ===========================================================================
% reverse gamma correction to R
R_double_norm = R_double./255;
if lt (R_double_norm, 0.04045)
    ori_R_double = ((R_double_norm)./(12.92));
else
    ori_R_double = ((R_double_norm + 0.055)./(1.055)).^2.4;
end
%ori_R_double = ori_R_double.*255;
ori_R_double = ori_R_double.*100;

% reverse gamma correction to G
G_double_norm = G_double./255;
if lt (G_double_norm, 0.04045)
    ori_G_double = ((G_double_norm)./(12.92));
else
    ori_G_double = ((G_double_norm + 0.055)./(1.055)).^2.4;
end
%ori_G_double = ori_G_double.*255;
ori_G_double = ori_G_double.*100;

% reverse gamma correction to B
B_double_norm = B_double./255;
if lt (B_double_norm, 0.04045)
    ori_B_double = ((B_double_norm)./(12.92));
else
    ori_B_double = ((B_double_norm + 0.055)./(1.055)).^2.4;
end
%ori_B_double = ori_B_double.*255;
ori_B_double = ori_B_double.*100;

% change gamma double RGB to uint8
ori_R_gamma = uint8(ori_R_double);
ori_G_gamma = uint8(ori_G_double);
ori_B_gamma = uint8(ori_B_double);

% new image in gamma sRGB (uint8)
gamma_RGB(:,:,1) = ori_R_gamma; % X plane
gamma_RGB(:,:,2) = ori_G_gamma; % Y plane
gamma_RGB(:,:,3) = ori_B_gamma; % Z plane

% ===========================================================================
% Conversion to XYZ color space
% ===========================================================================
% matrix alpha for conversion from sRGB to XYZ under D65 light source
alpha = [0.41239079926596 0.35758433938388 0.18048078840183; 
    0.21263900587151 0.71516867876776 0.07219231536073; 
    0.01933081871559 0.11919477979463 0.95053215224966];

% conversion from sRGB to XYZ for X, Y and Z plane
ori_X_double = ((ori_R_double)*alpha(1,1))+((ori_G_double)*alpha(1,2))+((ori_B_double)*alpha(1,3));
ori_Y_double = ((ori_R_double)*alpha(2,1))+((ori_G_double)*alpha(2,2))+((ori_B_double)*alpha(2,3));
ori_Z_double = ((ori_R_double)*alpha(3,1))+((ori_G_double)*alpha(3,2))+((ori_B_double)*alpha(3,3));

% ori_X_double, ori_Y_double and ori_Z_double are tristimulus values of the
% original image.

%% ORIGINAL IMAGE XYZ TO CIELAB

ref_X = 95.047;
ref_Y = 100.00;
ref_Z = 108.883;

var_ori_X = ori_X_double / ref_X;          %ref_X =  95.047   Observer= 2°, Illuminant= D65
var_ori_Y = ori_Y_double / ref_Y;          %ref_Y = 100.000
var_ori_Z = ori_Z_double / ref_Z;          %ref_Z = 108.883

for i = 1:1:N;
    for j = 1:1:M;
        if (var_ori_X(i,j) > 0.008856);
            var_ori_X(i,j) = var_ori_X(i,j)^(1/3);
        else
            var_ori_X(i,j) = ((1/3)*(29/6)*(29/6)*var_ori_X(i,j)) + (16/116);
        end
    end
end

for i = 1:1:N;
    for j = 1:1:M;
        if (var_ori_Y(i,j) > 0.008856); 
            var_ori_Y(i,j) = var_ori_Y(i,j)^(1/3);
        else
            var_ori_Y(i,j) = ((1/3)*(29/6)*(29/6)*var_ori_Y(i,j)) + (16/116);
        end
    end
end
    
for i = 1:1:N;
    for j = 1:1:M;
        if (var_ori_Z(i,j) > 0.008856);
            var_ori_Z(i,j) = var_ori_Z(i,j) ^ (1/3);
        else
            var_ori_Z(i,j) = ((1/3)*(29/6)*(29/6)*var_ori_Z(i,j)) + (16/116);
        end
    end
end
for i = 1:1:N;
    for j = 1:1:M;
        CIE_ori_L(i,j) = (116 * var_ori_Y(i,j)) - 16;
        CIE_ori_a(i,j) = 500 * (var_ori_X(i,j) - var_ori_Y(i,j));
        CIE_ori_b(i,j) = 200 * (var_ori_Y(i,j) - var_ori_Z(i,j));
    end
end

%% CHROMA, HUE AND INTENSITY

for i = 1:1:N;
    for j = 1:1:M;
            C(i,j) = (((CIE_ori_a(i,j))^2)+((CIE_ori_b(i,j))^2))^(1/2);  %Chroma
            H(i,j) = atan (abs((CIE_ori_b(i,j)/CIE_ori_a(i,j))));        %Hue
            H2(i,j) = H(i,j)*180*7/22;                                   %Change Hue from radians to degree
    end
end

old_L = uint8(CIE_ori_L);
old_a = uint8(CIE_ori_a);
old_b = uint8(CIE_ori_b);

% INTENSITY

for i = 1:1:N;
    for j = 1:1:M;
        I(i,j)= (1/3)*(R_double(i,j)+ G_double(i,j)+ B_double(i,j));
    end
end

%% CALCULATE NEW LUMINANCE
% 0 < alpha < 1
% alpha = (255 - intensity)/255
% alpha (i,j) = (255 - I(i,j))/255;
% L3(i,j) = CIE_ori_L(i,j)*(1+alpha(i,j)) - This is L2 in the paper;

for i = 1:1:N;
    for j = 1:1:M;
        alpha (i,j) = (255 - I(i,j))/255;
    end
end

for i = 1:1:N;
    for j = 1:1:M;
        L3(i,j) = CIE_ori_L(i,j)*(1 + alpha(i,j));
    end
end

new_L = uint8(L3);

%% RECOLORING PROCESS - TWO STAGES - LUMINANCE AND CHROMINANCE ENHANCEMENT
% Old luminance CIE_ori_L(i,j)
% New luminance L3(i,j)

for i = 1:1:N;
    for j = 1:1:M;
        C_ratio(i,j) = L3(i,j)/CIE_ori_L(i,j);
        C_new(i,j) = C_ratio(i,j)*C(i,j);   % New Chroma value
        if (CIE_ori_a(i,j) > 0);
            CIE_ori_a_bar(i,j) = C_new(i,j)*cos(H(i,j));  % New a value
        else
            CIE_ori_a_bar(i,j) = -C_new(i,j)*cos(H(i,j)); % New a value
        end
    end
end

for i = 1:1:N;
    for j = 1:1:M;
        if (CIE_ori_b(i,j) > 0);
            CIE_ori_b_bar(i,j) = C_new(i,j)*sin(H(i,j));  % New b value
        else
            CIE_ori_b_bar(i,j) = -C_new(i,j)*sin(H(i,j)); % New b value
        end
    end
end

new_a = uint8(CIE_ori_a_bar);
new_b = uint8(CIE_ori_b_bar);

%% LAB to XYZ

var_enh_Y = ((L3 + 16)/116);
var_enh_X = var_enh_Y +(CIE_ori_a_bar/500);
var_enh_Z = var_enh_Y -(CIE_ori_b_bar/200);

for i = 1:1:N;
    for j = 1:1:M;
       if ((var_enh_X(i,j)^3) > 0.008856);
            var_enh_X(i,j) = (var_enh_X(i,j))^3;
        else
            var_enh_X(i,j) = ((var_enh_X(i,j) - (16/116))/7.787);
        end
    end
end

for i = 1:1:N;
    for j = 1:1:M;
        if ((var_enh_Y(i,j)^3) > 0.008856);
            var_enh_Y(i,j) = var_enh_Y(i,j)^3;
        else
            var_enh_Y(i,j) = ((var_enh_Y(i,j) - (16/116))/7.787);
        end
    end
end

for i = 1:1:N;
    for j = 1:1:M;
        if ((var_enh_Z(i,j)^3) > 0.008856);
            var_enh_Z(i,j) = var_enh_Z(i,j)^3;
        else
            var_enh_Z(i,j) = ((var_enh_Z(i,j) - (16/116))/7.787);
        end
    end
end

%%
maoi_X_double = var_enh_X*ref_X;          %ref_X =  95.047   Observer= 2°, Illuminant= D65
maoi_Y_double = var_enh_Y*ref_Y;          %ref_Y = 100.000
maoi_Z_double = var_enh_Z*ref_Z;          %ref_Z = 108.883

% maoi_X_double, maoi_Y_double and maoi_Z_double are tristimulus values of the
% enhanced image.

% ===========================================================================
% Conversion to RGB color space
% ===========================================================================
% matrix inverse alpha for conversion from XYZ to sRGB under D65 light source
rev_alpha = [3.24096994190451,-1.53738317757009,-0.498610760292993;
    -0.969243636280870,1.87596750150771,0.0415550574071777;
    0.0556300796970030,-0.203976958888986,1.05697151424288];

% oonversion from XYZ to sRGB for R, G and B plane
maoigamma_R_double = ((maoi_X_double)*rev_alpha(1,1))+((maoi_Y_double)*rev_alpha(1,2))+((maoi_Z_double)*rev_alpha(1,3));
maoigamma_G_double = ((maoi_X_double)*rev_alpha(2,1))+((maoi_Y_double)*rev_alpha(2,2))+((maoi_Z_double)*rev_alpha(2,3));
maoigamma_B_double = ((maoi_X_double)*rev_alpha(3,1))+((maoi_Y_double)*rev_alpha(3,2))+((maoi_Z_double)*rev_alpha(3,3));

% ===========================================================================
% Gamma correction to the original RGB
% ===========================================================================
% gamma correction to R
maoinewgamma_R_double_norm = maoigamma_R_double./100;
if lt (maoinewgamma_R_double_norm, 0.0031308)
    maoi_R_double = 12.92*maoinewgamma_R_double_norm;
else
    maoi_R_double = ((1.055)*((maoinewgamma_R_double_norm).^(1/2.4)))-0.055;
end
maoi_R_double = maoi_R_double.*255;

% gamma correction to G
maoinewgamma_G_double_norm = maoigamma_G_double./100;
if lt (maoinewgamma_G_double_norm, 0.0031308)
    maoi_G_double = 12.92*maoinewgamma_G_double_norm;
else
    maoi_G_double = ((1.055)*((maoinewgamma_G_double_norm).^(1/2.4)))-0.055;
end
maoi_G_double = maoi_G_double.*255;

% gamma correction to B
maoinewgamma_B_double_norm = maoigamma_B_double./100;
if lt (maoinewgamma_B_double_norm, 0.0031308)
    maoi_B_double = 12.92*maoinewgamma_B_double_norm;
else
    maoi_B_double = ((1.055)*((maoinewgamma_B_double_norm).^(1/2.4)))-0.055;
end
maoi_B_double = maoi_B_double.*255;

% change new double RGB to uint8
maoi_R = uint8(maoi_R_double);
maoi_G = uint8(maoi_G_double);
maoi_B = uint8(maoi_B_double);

% remove any imaginary values in RGB
maoi_R = real(maoi_R);
maoi_G = real(maoi_G);
maoi_B = real(maoi_B);

% new image in sRGB (uint8)
maoi_RGB(:,:,1) = maoi_R; 
maoi_RGB(:,:,2) = maoi_G; 
maoi_RGB(:,:,3) = maoi_B; 

%%

figure (1)
subplot(1,2,1)
imshow (ori), title('Distorted Image');
subplot(1,2,2)
imshow (maoi_RGB), title('Enhanced Image');

imwrite(maoi_RGB,'improved.bmp');      %save enhanced image

