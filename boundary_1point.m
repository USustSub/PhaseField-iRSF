function [kglob,rhs] = boundary_1point(kglob,rhs,dof,displ)

% modify rhs
rhs(:,1) = rhs(:,1) - displ*kglob(:,dof);
rhs(dof,1) = displ;

% modify kglob
kglob(dof,:) = 0;
kglob(:,dof) = 0;
kglob(dof,dof) = 1;