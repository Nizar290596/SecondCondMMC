solve
(
     fvm::ddt(rho,z)
   + fvm::div(phi,z)
   - fvm::laplacian(rho * DEff,z)
);

z.max(0.0);
            
z.min(1.0);
