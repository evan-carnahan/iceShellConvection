function [D,G,I]=build_ops(Grid)
       % author: Evan Carnahan
       % date: 9/14/19
       % HW1 Q3 (2)
       % description:
       % This function computes the discrete divergence and gradient matrices on a
       % regular staggered grid using central difference approximations. The
       % discrete gradient assumes homogeneous boundary conditions.
       % Input:
       % Grid = structure containing all pertinent information about the grid.
       % Output:
       % D = Nx by Nx+1 discrete divergence matrix
       %old version
       %D = 1/Grid.dx * (spdiags(ones(Grid.dof_xmax,1),1,Grid.dof_xmax,Grid.dof_f_xmax)+ ...
       %    spdiags(-1*ones(Grid.dof_xmax,1),0,Grid.dof_xmax,Grid.dof_f_xmax));
       
       % G = Nx+1 by Nx discrete gradient matrix
       %G = -D';
       %G(Grid.dof_f_xmin,Grid.dof_xmin) = 0; G(Grid.dof_f_xmax,Grid.dof_xmax) = 0;
       % I = Nx by Nx identity matrix
       %I = spdiags(ones(Grid.dof_xmax,1),0,Grid.dof_xmax,Grid.dof_xmax);
       
       %%update 10/15 ->
       Dx1 = spdiags([-ones(Grid.Nx,1) ones(Grid.Nx,1)]/Grid.dx,[0 1],Grid.Nx,Grid.Nx+1);
       Iy = speye(Grid.Ny); 
       Dx = kron(Dx1,Iy);
       Ix = speye(Grid.Nx);
       if Grid.Ny ~= 1
           Dy1 = spdiags([-ones(Grid.Ny,1) ones(Grid.Ny,1)]/Grid.dy,[0 1],Grid.Ny,Grid.Ny+1);
           Dy = kron(Ix,Dy1);
           D = [Dx,Dy];
       else
           D = Dx;
       end
       I = kron(Ix,Iy);
       G = -D';
       
       %set x-fluxes to zero
       %Gx = G(1:Grid.Nfx,:);
       %Gy = G(Grid.Nfx+1:Grid.Nf,:);
       %G([Grid.dof_xmin,Grid.dof_xmax,Grid.dof_xmax,Grid.dof_xmin],...
       %    [Grid.dof_ymin,Grid.dof_ymin,Grid.dof_ymax,Grid.dof_ymax]) = 0;
       %set y-fluxes to zero
       
       %G(,Grid.dof_xmax) = 0;
       G(Grid.dof_f_xmin,Grid.dof_xmin) = 0; G(Grid.dof_f_xmax,Grid.dof_xmax) = 0;
       G(Grid.dof_f_ymin,Grid.dof_ymin) = 0; G(Grid.dof_f_ymax,Grid.dof_ymax) = 0;
       
       % Example call:
       % >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
       % >> Grid = build_grid(Grid);
       % >> [D,G,I]=build_ops(Grid);