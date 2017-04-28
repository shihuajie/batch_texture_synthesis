mex -Iinclude -IM:\local\boost_1_49_0 src/nnfMex.cpp -output dist/nnmex
mex -Iinclude -IM:\local\boost_1_49_0 src/voteMex.cpp -output dist/votemex
mex -Iinclude -IM:\local\boost_1_49_0 src/segmentmex.cpp -output dist/segmentmex
mex -Iinclude -IM:\local\boost_1_49_0 src/ggdtmex.cpp -output dist/ggdtmex
mex -Iinclude -IM:\local\boost_1_49_0 src/latticemex.cpp -output dist/latticemex
mex -Iinclude -IM:\local\boost_1_49_0 src/colorHistMex.cpp -output dist/colorhistmex
mex -Iinclude -IM:\local\boost_1_49_0 src/jigsawMex.cpp -output dist/jigsawMex

% debug mode
% mex -g -Iinclude -IC:\boost_1_49_0 src/nnfMex.cpp -output dist/nnmex
% mex -g -Iinclude -IC:\boost_1_49_0 src/voteMex.cpp -output dist/votemex
% mex -g -Iinclude -IC:\boost_1_49_0 src/segmentmex.cpp -output dist/segmentmex
% mex -g -Iinclude -IC:\boost_1_49_0 src/ggdtmex.cpp -output dist/ggdtmex
% mex -g -Iinclude -IC:\boost_1_49_0 src/latticemex.cpp -output dist/latticemex 
% mex -g -Iinclude -IC:\boost_1_49_0 src/colorHistMex.cpp -output dist/colorhistmex
% mex -g -Iinclude -IC:\boost_1_49_0 src/jigsawMex.cpp -output dist/jigsawMex