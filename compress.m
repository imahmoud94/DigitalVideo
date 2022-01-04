function cframe = compress(iframe, dct)

    c = iframe;
    imgycbcr = rgb2ycbcr(iframe);
    Y = imgycbcr(:,:,1);
    Cb = imgycbcr(:,:,2);
    Cr = imgycbcr(:,:,3);

    
    D = dct;
    T = dctmtx(D);

    Y = im2double(Y);
    Cb = im2double(Cb);
    Cr = im2double(Cr);
    BY = blkproc(Y,[D D],'P1*x*P2',T,T');
    BCb = blkproc(Cb,[D D],'P1*x*P2',T,T');
    BCr = blkproc(Cr,[D D],'P1*x*P2',T,T');

    mask4 = [1 1 0 0
            1 0 0 0
            0 0 0 0
            0 0 0 0];

    mask8 = [1   1   1   1   0   0   0   0
            1   1   1   0   0   0   0   0
            1   1   0   0   0   0   0   0
            1   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0];

    mask16 = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
            1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
            1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
            1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
            1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0
            1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
            1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    if(dct == 4)
		B2Y = blkproc(BY,[D D],'P1.*x',mask4);
		B2Cb = blkproc(BCb,[D D],'P1.*x',mask4);
		B2Cr = blkproc(BCr,[D D],'P1.*x',mask4);
	elseif(dct == 8)
		B2Y = blkproc(BY,[D D],'P1.*x',mask8);
		B2Cb = blkproc(BCb,[D D],'P1.*x',mask8);
		B2Cr = blkproc(BCr,[D D],'P1.*x',mask8);
	else
		B2Y = blkproc(BY,[D D],'P1.*x',mask16);
		B2Cb = blkproc(BCb,[D D],'P1.*x',mask16);
		B2Cr = blkproc(BCr,[D D],'P1.*x',mask16);
	end
    
    I2Y = blkproc(B2Y,[D D],'P1*x*P2',T',T);
    I2Cb = blkproc(B2Cb,[D D],'P1*x*P2',T',T);
    I2Cr = blkproc(B2Cr,[D D],'P1*x*P2',T',T);

    imgycbcr(:,:,1) = im2uint8(I2Y);
    imgycbcr(:,:,2) = im2uint8(I2Cb);
    imgycbcr(:,:,3) = im2uint8(I2Cr);
    cframe = ycbcr2rgb(imgycbcr);

end