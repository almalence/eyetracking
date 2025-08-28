function [I, mapOrAlpha] = buffered2im( jImage, varargin )
%BUFFERED2IM Convert a Java BufferedImage to a Matlab image
%
%    I = BUFFERED2IM(JIMAGE) converts the Java BufferedImage, JIMAGE to the
%    Matlab image, I
%
%    I = BUFFERED2IM(JIMAGE,BGCOLOR) converts the Java BufferedImage,
%    JIMAGE whose background color is BGCOLOR, [red,green,blue] (double, 0 to 1)
%    to the Matlab image, I.  
%
%   [I,MAP] = BUFFERED2IM(JIMAGE) additionaly returns the colormap, MAP if
%   JIMAGE is of type BufferedImage.TYPE_BYTE_INDEXED.
%
%   [I,ALPHA] = BUFFERED2IM(...) additionaly returns the alpha channel,
%   ALPHA 
%
%   Example
%   -------
%   This example reads an image into the MATLAB workspace and then uses
%   im2buffered to convert it into an instance of the Java BufferedImage 
%   class and than reconverts the Java BufferedImage to a Matlab image
%
%   I = imread('moon.tif');
%   [jImage,g2d] = im2Buffered(I);
%   ...manipulate image using g2d...
%   I = buffered2im(jImage);
%
%   Input-output specs
%   ------------------ 
%   JIMAGE:   java.awt.BufferedImage object
%
%   BGCOLOR:  RGB triplet, double between 0 and 1
%
%   I:        2-D, or 3-D, real, full matrix
%             uint8, uint16, or double
%
%   MAP:      2-D, real, full matrix
%             size(MAP,1) <= 256
%             size(MAP,2) == 3
%             double 
%
%  ALPHA:     2-D, real, full matrix
%             same dimension as I
%             uint32
%
%  Code for converting Java image data buffer to RGB triplets is taken from
%  Altman, Yair. 'Undocumented Secrets of Matlab-Java Programming'
%
%  Collin Pecora 7/2018
%
%  See also im2buffered, im2java

    % Don't run on platforms with incomplete Java support
    error(javachk('awt','buffered2im'));
    % Check number of inputs and outputs
    narginchk(1,2);
    nargoutchk(0,2);
    % Pre-allocate
    mapOrAlpha = [];

    if isa(jImage,'java.awt.image.BufferedImage')
        % Image width
        imgWidth = jImage.getWidth;
        % Image height
        imgHeight = jImage.getHeight;  
        % Color model
        jColorModel = jImage.getColorModel;
        % Flag indicating presence of alpha channel, pixels may still be
        % all opaque, checking for an all opaque alpha channel is slower
        % than allowing applyBackground() to execute
        hasAlpha = jColorModel.hasAlpha();
        % Check for optional background color input argument
        if nargin == 2
            bgColor = varargin{1};
            % Validate background color
            validateattributes(bgColor,...
                {'double'},...
                {'row','size',[1,3],'nonempty','>=',0,'<=',1},...
                'buffered2im')
            % Convert background color to uint8 and java.awt.Color
            bgColor = uint8(round(bgColor*255)); 
            bgColor = java.awt.Color(bgColor(1),bgColor(2),bgColor(3),uint8(255));
            if hasAlpha
                applyBackground()
            end
        end        
        
        switch jImage.getType
            
            case java.awt.image.BufferedImage.TYPE_INT_ARGB 
            % (2) Represents an image with 8-bit RGBA color components
            % packed into integer pixels. The image has a DirectColorModel
            % with alpha. The color data in this image is considered not to
            % be premultiplied with alpha. When this type is used as the
            % imageType argument to a BufferedImage constructor, the
            % created image is consistent with images created in the JDK1.1
            % and earlier releases.
            
                handleDefaultImage();
                
%                 pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
%                 
%                 if hasAlpha
%                     pixels = reshape(pixels,[4,imgWidth*imgHeight]);
%                     mapOrAlpha = reshape(pixels(4,:),[imgWidth,imgHeight])';
%                 else
%                     pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
%                     pixels = reshape(pixels,[3,imgWidth*imgHeight]);
%                 end
%                 
%                 R = reshape(pixels(3,:),[imgWidth,imgHeight])';
%                 B = reshape(pixels(2,:),[imgWidth,imgHeight])';
%                 G = reshape(pixels(1,:),[imgWidth,imgHeight])'; 
% 
%                 I = cat(3,R,B,G);  
                
            case java.awt.image.BufferedImage.TYPE_4BYTE_ABGR
            % (6) Represents an image with 8-bit RGBA color components with
            % the colors Blue, Green, and Red stored in 3 bytes and 1 byte
            % of alpha. The image has a ComponentColorModel with alpha.
            % The color data in this image is considered not to be
            % premultiplied with alpha. The byte data is interleaved in a
            % single byte array in the order A, B, G, R from lower to
            % higher byte addresses within each pixel    
 
                pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
                pixels = reshape(pixels,[4,imgWidth*imgHeight]);
                
                R = reshape(pixels(4,:),[imgWidth,imgHeight])';
                B = reshape(pixels(3,:),[imgWidth,imgHeight])';
                G = reshape(pixels(2,:),[imgWidth,imgHeight])';
                
                mapOrAlpha = reshape(pixels(1,:),[imgWidth,imgHeight])';

                I = cat(3,R,B,G); 
                
            case java.awt.image.BufferedImage.TYPE_3BYTE_BGR
            % (5) Represents an image with 8-bit RGB color components,
            % corresponding to a Windows-style BGR color model) with the
            % colors Blue, Green, and Red stored in 3 bytes. There is no
            % alpha. The image has a ComponentColorModel. When data with
            % non-opaque alpha is stored in an image of this type, the
            % color data must be adjusted to a non-premultiplied form and
            % the alpha discarded, as described in the 
            % java.awt.AlphaComposite documentation. 

                pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
                pixels = reshape(pixels,[3,imgWidth*imgHeight]);

                R = reshape(pixels(3,:),[imgWidth,imgHeight])';
                B = reshape(pixels(2,:),[imgWidth,imgHeight])';
                G = reshape(pixels(1,:),[imgWidth,imgHeight])';
                
                I = cat(3,R,B,G); 
                
            case java.awt.image.BufferedImage.TYPE_INT_RGB

                pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
                pixels = reshape(pixels,[4,imgWidth*imgHeight]);

                R = reshape(pixels(3,:),[imgWidth,imgHeight])';
                B = reshape(pixels(2,:),[imgWidth,imgHeight])';
                G = reshape(pixels(1,:),[imgWidth,imgHeight])';
                
                I = cat(3,R,B,G);                 
                
            case java.awt.image.BufferedImage.TYPE_BYTE_INDEXED
            % (13) Represents an indexed byte image. When this type is used
            % as the imageType argument to the BufferedImage constructor
            % that takes an imageType argument but no ColorModel argument,
            % an IndexColorModel is created with a 256-color 6/6/6 color
            % cube palette with the rest of the colors from 216-255
            % populated by grayscale values in the default sRGB ColorSpace.
            % When color data is stored in an image of this type, the
            % closest color in the colormap is determined by the
            % IndexColorModel and the resulting index is stored.
            % Approximation and loss of alpha or color components can result,
            % depending on the colors in the IndexColorModel colormap.
                
                try
                    % Map
                    colorCount = jColorModel.getMapSize();
                    
%                     sz = jColorModel.getComponentSize;

                    getReds = io.JavaMethodWrapper('java.awt.image.IndexColorModel','getReds(byte[])');
                    getGreens = io.JavaMethodWrapper('java.awt.image.IndexColorModel','getGreens(byte[])');
                    getBlues = io.JavaMethodWrapper('java.awt.image.IndexColorModel','getBlues(byte[])');                    

                    reds = zeros(1,colorCount,'int8');
                    greens = zeros(1,colorCount,'int8');
                    blues = zeros(1,colorCount,'int8');                   

                    [~,reds] = getReds.invoke(jColorModel,reds);
                    [~,greens] = getGreens.invoke(jColorModel,greens);
                    [~,blues] = getBlues.invoke(jColorModel,blues);                    

                    mapOrAlpha = double([...
                        typecast(reds,'uint8'),...
                        typecast(greens,'uint8'),...
                        typecast(blues,'uint8')])./255;
                    
                    if hasAlpha
                        getAlphas = io.JavaMethodWrapper('java.awt.image.IndexColorModel','getAlphas(byte[])');
                        alphas = zeros(1,colorCount,'int8');
                        [~,alphas] = getAlphas.invoke(jColorModel,alphas);
                        alphas = double(typecast(alphas,'uint8'))./255;
                        
                        mapOrAlpha = [mapOrAlpha,alphas];
                    end
                catch
                end
                
                
                
                I = typecast(jImage.getRaster.getDataBuffer.getData,'uint8');
                I = reshape(I,[imgWidth,imgHeight])';    
                
            case java.awt.image.BufferedImage.TYPE_BYTE_GRAY                
            % (10) Represents a unsigned byte grayscale image, non-indexed.
            % This image has a ComponentColorModel with a CS_GRAY
            % java.awt.color.ColorSpace. When data with non-opaque alpha is
            % stored in an image of this type, the color data must be
            % adjusted to a non-premultiplied form and the alpha discarded,
            % as described in the java.awt.AlphaComposite documentation.
            
                I = typecast(jImage.getRaster.getDataBuffer.getData,'uint8');
                I = reshape(I,[imgWidth,imgHeight])';  
                
            case java.awt.image.BufferedImage.TYPE_USHORT_GRAY
            % (11)  Represents an unsigned short grayscale image, non-indexed).
            % This image has a ComponentColorModel with a CS_GRAY ColorSpace.
            % When data with non-opaque alpha is stored in an image of this
            % type, the color data must be adjusted to a non-premultiplied
            % form and the alpha discarded, as described in the
            % java.awt.AlphaComposite documentation. 
            
                I = typecast(jImage.getRaster.getDataBuffer.getData,'uint16');
                I = reshape(I,[imgWidth,imgHeight])';
                
            case java.awt.image.BufferedImage.TYPE_CUSTOM
            % Unkown or custom type. Convert to a TYPE_INT_ARGB 
            
                jNewImage = java.awt.image.BufferedImage(jImage.getWidth,...
                    jImage.getHeight,java.awt.image.BufferedImage.TYPE_INT_ARGB);

                g2d = jNewImage.createGraphics();
                g2d.drawImage(jImage, 0, 0, jImage.getWidth(), jImage.getHeight(), []);
                g2d.dispose(); 

                jImage = jNewImage;

                handleDefaultImage();
                
            otherwise
                error('buffered2im:WrongType','Unsupported BufferedImage type') 
        end
        
    else
        error('buffered2im:InvalidInput',...
            'First input argument must be a java.awt.image.BufferedImage')        
    end
    
    function applyBackground()
    % APPLYBACKGROUND draws jImage onto a temporary buffered image whose
    % background color is bgColor and than sets jImage to the temporary image
    % The defualt AlphaCompositing rule is applied
    % (AlphaComposite.SRC_OVER)
    % See java.awt.AlphaComposite
    
        % Create a temporary BufferedImage, same size and type as jImage
        temp = java.awt.image.BufferedImage(imgWidth,imgHeight,jImage.getType);
        % Get Graphics2D object
        g2d = temp.createGraphics();
        % Set fill color to bgColor
        g2d.setPaint(bgColor);
        % Fill background with bgColor
        g2d.fillRect(0,0,imgWidth,imgHeight);
        % Draw jImage onto temp
        g2d.drawImage(jImage,0,0,[]);
        % Dispose graphic object
        g2d.dispose();
        jImage = temp;
    end

    function handleDefaultImage()
        
        pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';

        if hasAlpha
            pixels = reshape(pixels,[4,imgWidth*imgHeight]);
            mapOrAlpha = reshape(pixels(4,:),[imgWidth,imgHeight])';
        else
            pixels = typecast(jImage.getRaster.getDataBuffer.getData,'uint8')';
            pixels = reshape(pixels,[3,imgWidth*imgHeight]);
        end

        R = reshape(pixels(3,:),[imgWidth,imgHeight])';
        B = reshape(pixels(2,:),[imgWidth,imgHeight])';
        G = reshape(pixels(1,:),[imgWidth,imgHeight])'; 

        I = cat(3,R,B,G); 
    end
end
% // Image type is not recognized so it must be a customized image. This type is only used as a return value for the getType() method.
% public static final int TYPE_CUSTOM = 0;
% 
% // Represents an image with 8-bit RGB color components packed into integer pixels. The image has a DirectColorModel without alpha. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_INT_RGB = 1;
% 
% // Represents an image with 8-bit RGBA color components packed into integer pixels. The image has a DirectColorModel with alpha. The color data in this image is considered not to be premultiplied with alpha. When this type is used as the imageType argument to a BufferedImage constructor, the created image is consistent with images created in the JDK1.1 and earlier releases.
% public static final int TYPE_INT_ARGB = 2;
% 
% // Represents an image with 8-bit RGBA color components packed into integer pixels. The image has a DirectColorModel with alpha. The color data in this image is considered to be premultiplied with alpha.
% public static final int TYPE_INT_ARGB_PRE = 3;
% 
% // Represents an image with 8-bit RGB color components, corresponding to a Windows- or Solaris- style BGR color model, with the colors Blue, Green, and Red packed into integer pixels. There is no alpha. The image has a DirectColorModel. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_INT_BGR = 4;
% 
% // Represents an image with 8-bit RGB color components, corresponding to a Windows-style BGR color model) with the colors Blue, Green, and Red stored in 3 bytes. There is no alpha. The image has a ComponentColorModel. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_3BYTE_BGR = 5;
% 
% // Represents an image with 8-bit RGBA color components with the colors Blue, Green, and Red stored in 3 bytes and 1 byte of alpha. The image has a ComponentColorModel with alpha. The color data in this image is considered not to be premultiplied with alpha. The byte data is interleaved in a single byte array in the order A, B, G, R from lower to higher byte addresses within each pixel.
% public static final int TYPE_4BYTE_ABGR = 6;
% 
% // Represents an image with 8-bit RGBA color components with the colors Blue, Green, and Red stored in 3 bytes and 1 byte of alpha. The image has a ComponentColorModel with alpha. The color data in this image is considered to be premultiplied with alpha. The byte data is interleaved in a single byte array in the order A, B, G, R from lower to higher byte addresses within each pixel.
% public static final int TYPE_4BYTE_ABGR_PRE = 7;
% 
% // Represents an image with 5-6-5 RGB color components (5-bits red, 6-bits green, 5-bits blue) with no alpha. This image has a DirectColorModel. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_USHORT_565_RGB = 8;
% 
% // Represents an image with 5-5-5 RGB color components (5-bits red, 5-bits green, 5-bits blue) with no alpha. This image has a DirectColorModel. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_USHORT_555_RGB = 9;
% 
% // Represents a unsigned byte grayscale image, non-indexed. This image has a ComponentColorModel with a CS_GRAY java.awt.color.ColorSpace. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_BYTE_GRAY = 10;
% 
% // Represents an unsigned short grayscale image, non-indexed). This image has a ComponentColorModel with a CS_GRAY ColorSpace. When data with non-opaque alpha is stored in an image of this type, the color data must be adjusted to a non-premultiplied form and the alpha discarded, as described in the java.awt.AlphaComposite documentation.
% public static final int TYPE_USHORT_GRAY = 11;
% 
% // Represents an opaque byte-packed 1, 2, or 4 bit image. The image has an IndexColorModel without alpha. When this type is used as the imageType argument to the BufferedImage constructor that takes an imageType argument but no ColorModel argument, a 1-bit image is created with an IndexColorModel with two colors in the default sRGB ColorSpace: {0, 0, 0} and {255, 255, 255}. Images with 2 or 4 bits per pixel may be constructed via the BufferedImage constructor that takes a ColorModel argument by supplying a ColorModel with an appropriate map size. Images with 8 bits per pixel should use the image types TYPE_BYTE_INDEXED or TYPE_BYTE_GRAY depending on their ColorModel. When color data is stored in an image of this type, the closest color in the colormap is determined by the IndexColorModel and the resulting index is stored. Approximation and loss of alpha or color components can result, depending on the colors in the IndexColorModel colormap.
% public static final int TYPE_BYTE_BINARY = 12;
% 
% // Represents an indexed byte image. When this type is used as the imageType argument to the BufferedImage constructor that takes an imageType argument but no ColorModel argument, an IndexColorModel is created with a 256-color 6/6/6 color cube palette with the rest of the colors from 216-255 populated by grayscale values in the default sRGB ColorSpace. When color data is stored in an image of this type, the closest color in the colormap is determined by the IndexColorModel and the resulting index is stored. Approximation and loss of alpha or color components can result, depending on the colors in the IndexColorModel colormap.
% public static final int TYPE_BYTE_INDEXED = 13;

