function f_releaseCamera(vid)
%%% Close preview if open:
if exist('vid','var') == 1
  if ~isempty(vid) && strcmp(vid.previewing,'on') % vid not empty means 
    % there is an image source
    closepreview(vid); % Gets closed in case it was already opened
  end
  delete(vid); % Delete vid
end

end

