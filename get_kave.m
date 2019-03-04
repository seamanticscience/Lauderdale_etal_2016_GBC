function kave=get_kave(depth_axis,mixedlayer_depth,varargin)
%function kave=get_kave(depth_axis,mixedlayer_depth)
%
% Returns a 2d grid of z-levels associated with the depth of the mixed
% layer. If mixedlayer_depth is a single, arbitrary depth value, then kave
% will return a 2d array of constant values.
%
% JML Dec 2013
global grid

if nargin>2
    nx=varargin{1};
    ny=varargin{2};
    nt=size(mixedlayer_depth,3);
else
    nx=grid.nx;
    ny=grid.ny;
    nt=size(mixedlayer_depth,3);
end

if length(mixedlayer_depth(:))==1
    kave=ones(nx,ny,nt).*find(depth_axis>=mixedlayer_depth,1,'first');
else
    kave=ones(nx,ny,nt);
    for t=1:nt;
        for i=1:nx
            for j=1:ny
                if ~isnan(mixedlayer_depth(i,j,t)) && mixedlayer_depth(i,j,t)>=depth_axis(1)
                    kave(i,j,t)=find(depth_axis>=mixedlayer_depth(i,j,t),1,'first');
                else
                    kave(i,j,t)=1;
                end
            end
        end
    end
end

end