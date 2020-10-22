%function HaarTransformationMatrix=ConstructHaarWaveletTransformationMatrix(WidthOfSquareMatrix)
%--------------------------------------------------------------------------
%Create Haarwavelet transformation matrix H for the matrix vector
%mulplication implimentation of Haar wavelet transformation.
%This function uses the following nice formula to create the Haar
%transformation matrix:
%               H_n=1/sqrt(2)[H_(n/2) kron (1 1)
%                             I_(n/2) kron (1 -1)],
%                              where 'kron' denotes the kronecker product.
%The iteration starts with H_1=[1]. The normalization constant 1/sqrt(2)
%ensure that H_n^T*H_n=I, where I is identity matrix. Haar wavelets are the
%rows of H_n.
%--------------------------------------------------------------------------
%function HaarTransformationMatrix=ConstructHaarWaveletTransformationMatrix(WidthOfSquareMatrix)
% Input:
%       WidthOfSquareMatrix: the width of sqaure Haar wavelet
%                                  transformation matrix. 
% Output:
%       HaarTransformationMatrix: Ceated Haar transformation matrix and
%       it's size is the power of 2,i.e., 2, 4,
%       8,16,32,64,128,256,512,1024,2048,4096,etc.

function [h] = haar(n)
h = [1];
if n > 2
    h = haar(n/2);
end
% calculate upper haar part
h_n = kron(h,[1,1]); 
% calculate lower haar part 
h_i = kron(eye(length(h)),[1,-1]);
% combine parts
h = [h_n; h_i]/sqrt(2);
end

%disp(haar(8));