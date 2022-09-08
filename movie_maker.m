%% Image Writer
clear all;
close all;
clc;
% load the images
num_images = 1101;
 images    = cell(num_images,1);
 for i=1:num_images
     images{i} = imread(sprintf('Img/img_%d.png',i));
 end
 % create the video writer with 30 fps
 writerObj = VideoWriter('Velocity.avi');
 writerObj.FrameRate = 30;
   % open the video writer
   open(writerObj);
   % write the frames to the video
    for u=1:num_images
       % convert the image to a frame
       frame = im2frame(images{u});
           writeVideo(writerObj, frame);

   end
   % close the writer object
   close(writerObj);
   implay('Velocity.avi');