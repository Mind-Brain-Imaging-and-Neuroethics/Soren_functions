function [output] = t3D(arrayin)

output = arrayin;

for c = 1:size(arrayin,3)
    output(:,:,c) = squeeze(arrayin(:,:,c))';
end