% Function to perform image stitching in a fully automated manner
function [outputs] = StitchAuto(sInfoL, sInfoR)


corners = STI_FUNC.GetCorners(sInfoL, sInfoR);


outputs = [];