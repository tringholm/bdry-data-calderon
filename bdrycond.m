function bcMatrix = bdrycond(location,state)
% bcMatrix = 1-location.x;
bcMatrix = max(0,location.x);
% if x > 0
%     bcMatrix = 1-location.x - location.y;
% end
end