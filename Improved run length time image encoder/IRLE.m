clear all
    clc
    %%%%%%%                Encoding             %%%%%%%

    % % Reading image such that the image name and extension must be filled
    % % as 'image_name.extension' , and the image must be in project folder.

    img=imread('1.png');

    % % Detemining image dimension

    [d1,d2,d3] = size(img);

    % % Data type conversion.

    M=double(img(:));

    % % Initializing vectors.

    [P,P_l,occ,reimg]=zeros;

    % % Finding the Difference vactor.

    for i= 2:length(M)
        P(i-1)= M(i)-M(i-1);
    end

    % % Converting P into logical format such that
    % % non zero elements are denoted by 1 , and the
    % % zero elements are denoted by 0.

    for i=1:length(P)
        if P(i)~=0
            P_l(i)=1;
        end
        if i == length(P)
            P_l(i+1)=1;
        end
    end

    % % Indexing non zero elements

    index=find(P_l);

    % % Finding the corresponding values in M

    values=M(index);

    % % Constructing the occurance vector

    occ(1)=index(1);
    for i=2:length(index)
        occ(i)=index(i)-index(i-1);
    end
    occ=occ(:);

    % % Concatenating the values into 1 matrix

    conc=[values,occ];

    % % Saving compressed image

    save('Compressed image','conc')

    % % Compression ratio by comparing the length of both 
    % %  the reduced image values array and the original image array.

    CompressionRatio=(1-(length(values)/length(M)))*100

    %%%%%%%                  Decoding             %%%%%%%

    % %  Reconstructing image array from concatenated pairs
    for i=1:length(values)
        if i==1
            j=0;
        else
            j=length(reimg);
        end
        for k=1:occ(i)
            reimg(j+k)=values(i);
        end
    end

    reimg=reimg(:);
    % % Converting data types from double to 8 bit integer values
    % % and constructing image matrix from the previously constructed 
    % % image array.

    const_img=uint8(reshape(reimg,d1,d2,d3));

    % %  Previewing both the original image and the reconstructed image

    figure
    subplot(1,2,1)
    imshow(img)
    subplot(1,2,2)
    imshow(const_img)
