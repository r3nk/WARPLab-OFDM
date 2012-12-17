%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OFDM Framework
%
% Copyright (c) 2012 Robin Klose 
%
% This script stops all transmissions on a WARP board running WARPLab 4x4. 
% It may be used to disable continuous transmissions or to perform a soft
% reset if something goes wrong. 
%
% The following parameters should be given: 
%
% node_id:     [1, 2] The node that will be stopped. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function warplab_4x4_stop(node_id)

% Load some global definitions:
warplab_defines

% Create socket handles and initialize nodes:
[socketHandles, packetNum] = warplab_initialize;

% Socket handle for the magic SYNC: 
udp_sync = socketHandles(1);

% Open two handles for two connected boards:
udp_node1 = socketHandles(2);
udp_node2 = socketHandles(3);

if node_id == 1
    udp_node = udp_node1;
elseif node_id == 2
    udp_node = udp_node2;
else
    error ('Bad node_id parameter: %d', node_id);
end

% Disable radio resources:
warplab_sendCmd(udp_node, TX_STOP, packetNum); 
warplab_sendCmd(udp_node, RADIO1TXBUFF_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO1_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO2TXBUFF_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO2_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO3TXBUFF_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO3_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO4TXBUFF_TXDIS, packetNum);
warplab_sendCmd(udp_node, RADIO4_TXDIS, packetNum);
