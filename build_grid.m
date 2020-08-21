function [Grid] = build_grid(Grid)
% author: Evan Carnahan
% Description:
% This function computes takes in minimal definition of the computational
% domain and grid and computes all pertinent information
% about the grid.
%
% Input:
% Grid.xmin = left boundary of the domain
% Grid.xmax = right bondary of the domain
% Grid.Nx   = number of grid cells

% Grid.xf = Nfx by 1 column vector of cell face locations
Grid.xf = (linspace(Grid.xmin,Grid.xmax,Grid.Nx+1))';
% Grid.dx = scalar cell width
Grid.dx = Grid.xf(2)-Grid.xf(1);
% Grid.xc = Nx by 1 column vector of cell center locations
Grid.xc = (Grid.xf(1:end-1)+Grid.dx/2); 
% Grid.Lx = scalar length of the domain
Grid.Lx = Grid.xmax-Grid.xmin;

%info for y-corrdinate 
if ~isfield(Grid,'ymin'); Grid.ymin = 0; end
if ~isfield(Grid,'ymax'); Grid.ymax = 1; end
if ~isfield(Grid,'Ny'); Grid.Ny = 1; end
Grid.yf = (linspace(Grid.ymin,Grid.ymax,Grid.Ny+1))';
Grid.dy = Grid.yf(2)-Grid.yf(1); 
Grid.Ly = Grid.ymax-Grid.ymin;
Grid.yc = (Grid.yf(1:end-1)+Grid.dy/2);

%update for 2D operations
Grid.N = Grid.Nx*Grid.Ny;
% Grid.Nfx = number of fluxes in x-direction
Grid.Nfx = Grid.N+Grid.Ny;
Grid.Nfy = Grid.N+Grid.Nx;
Grid.Nf = Grid.Nfx + Grid.Nfy;

% Grid.dof = Nx by 1 column vector from 1 to N containig the degrees of freedom, i.e. cell number
Grid.dof_f_x = reshape(1:Grid.Nfx,Grid.Ny,Grid.Nx+1);
Grid.dof_f_y = reshape(1:Grid.Nfy,Grid.Ny+1,Grid.Nx);
Grid.dof = reshape(1:Grid.N,Grid.Ny,Grid.Nx);

% Grid.dof_xmin  = scalar cell degree of freedom corrsponding to the left boundary
% Grid.dof_xmax  = scalar cell degree of freedom corrsponding to the right boundary
% Grid.dof_f_xmin = scalar face degree of freedom corrsponding to the left boundary
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right boundary
Grid.dof_f_xmin = Grid.dof_f_x(:,1);
Grid.dof_f_xmax = Grid.dof_f_x(:,end);
Grid.dof_f_ymin = Grid.Nfx+Grid.dof_f_y(1,:)';
Grid.dof_f_ymax = Grid.Nfx+Grid.dof_f_y(end,:)';

Grid.dof_xmin = Grid.dof(:,1);
Grid.dof_xmax = Grid.dof(:,end);
Grid.dof_ymin = Grid.dof(1,:)';
Grid.dof_ymax = Grid.dof(end,:)';


%1D area and volume, assumes square grid
if ~isfield(Grid,'zmin'); Grid.zmin = 0; end
if ~isfield(Grid,'zmax'); Grid.zmax = 1; end
if ~isfield(Grid,'Nz'); Grid.Nz = 1; end
Grid.zf = (linspace(Grid.zmin,Grid.zmax,Grid.Nz+1))';
Grid.dz = Grid.zf(2)-Grid.zf(1); 
Grid.Lz = Grid.zmax-Grid.zmin;
Grid.zc = (Grid.zf(1:end-1)+Grid.dz/2);

Grid.A = [ones(1,Grid.Nfx)*Grid.dy*Grid.dz ones(1,Grid.Nfy)*Grid.dx*Grid.dz]';
Grid.V = ones(Grid.N,1)*Grid.dx*Grid.dy*Grid.dz;


if ~isfield(Grid,'psi_x0');  Grid.psi_x0 = 'xmin_ymin'; end
if ~isfield(Grid,'psi_dir'); Grid.psi_dir = 'xy'; end
