function makeGif(F, nameString) 
  filename = nameString;
  for n = 1:size(F,2)
      im = frame2im(F(n)); 
      [imind,cm] = rgb2ind(im,300); 

      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.01); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.01); 
      end 
  end
end

