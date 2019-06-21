----------------------------------util.ns.CreateApproxSpace----------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	
--
--   Author: Josef Dubsky, Andreas Vogel
--
--------------------------------------------------------------------------------

PluginRequired("NavierStokes")

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("desilination_util.lua")
ug_load_script("util/conv_rates_static.lua")

numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 2, "number of prerefinements (parallel)")
bConvRates  = util.HasParamOption("-convRate", "compute convergence rates")
bBenchmarkRates = util.HasParamOption("-benchRate", "compute benchmark rates")

bStokes 	= util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-nolaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactjac", "If defined, exact jacobian FctCmpused")
bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type")
stab        = util.GetParam("-stab", "flow", "Stabilization type")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type")

porder = 1
vorder = 1
corder = 1
discType = "fe" 

local Viscosity	= 1e-3
local Um = 0.3 --max Inlet velocity
local H 
local B 

gridName = "cylinder"--"simpleStokes1806" 

local walls

if gridName == "wholeMesh" or gridName == "finalMesh" or gridName == "shiftMesh" then
	walls = "wall"
	H = 5
	B = 10
	discType = "fe"
	Um = 10
	dim = 3
elseif gridName == "mesh1104" then
	walls = "cylinder1,cylinder2,wall" 
	H = 2.828427125
	B = 1.25
	discType = "fe"
	Um = 10
	dim = 3
elseif gridName == "quadrefi" or gridName == "quadrefiPre" then
	H = 0.2
	B = H
	walls = "UpperWall,LowerWall,Cylinder1,Cylinder2,Cylinder3"
	dim = 2
	discType = "fe"
elseif gridName == "simpleStokes1806" then
	H = 0.5
	walls = "Wall"
	dim = 2
	discType = "fv"
elseif gridName == "cylinder" then
	H = 0.41
	walls = "UpperWall, LowerWall, CylinderWall"
	dim = 2
	discType = "fe"
else
	print("no such geometry")
end


-- Lets write some info about the choosen parameter
print(" Choosen Parater:")
print("    dim              = " .. dim)
print("    numTotalRefs     = " .. numRefs)
print("    numPreRefs       = " .. numPreRefs)
print("    grid             = " .. gridName)
print("    porder           = " .. porder)
print("    vorder           = " .. vorder)
print("    type             = " .. discType)
print("    only stokes      = " .. tostring(bStokes))
print("    no laplace       = " .. tostring(bNoLaplace))
print("    exact jacobian   = " .. tostring(bExactJac))
print("    peclet blend     = " .. tostring(bPecletBlend))
print("    upwind           = " .. upwind)
print("    stab             = " .. stab)
print("    diffLength       = " .. diffLength)

--------------------------------------------------------------------------------
-- Loading Domain and Domain Refinement
--------------------------------------------------------------------------------

function CreateDomain()

	InitUG(dim, AlgebraType("CPU", 1))
	local dom = Domain()
	LoadDomain(dom, "grids/"..gridName..".ugx")
	
	

	-- Create a refiner instance. This is a factory method
	-- which automatically creates a parallel refiner if required.
	local refiner =  GlobalDomainRefiner(dom)
	
	write("Pre-Refining("..numPreRefs.."): ")
	for i=1,numPreRefs do write(i .. " ");	refiner:refine(); end
	write("done. Distributing...")
	if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
		print("Error while Distributing Grid. Aborting.")
		exit();
	end
	write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
	for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end
	write("done.\n")
	
	--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
	
	return dom
end

function CreateApproxSpace(dom, discType, vorder, porder,corder)

	local approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder,corder)
	
	-- print statistic on the distributed dofs
	--approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
	--approxSpace:print_local_dof_statistic(2)
	
	return approxSpace
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------
globalNSDisc = nil
function CreateDomainDisc(approxSpace, discType, vorder, porder)
	
	local FctCmp -- = approxSpace:names()
	if 		dim == 1 then FctCmp = {"u", "p"}
	elseif  dim == 2 then FctCmp = {"u", "v", "p"}
	elseif  dim == 3 then FctCmp = {"u", "v", "w", "p"}
	else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end

	NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
	NavierStokesDisc:set_exact_jacobian(bExactJac)
	NavierStokesDisc:set_stokes(bStokes)
	NavierStokesDisc:set_laplace( not(bNoLaplace) )
	NavierStokesDisc:set_kinematic_viscosity( Viscosity );
	globalNSDisc = NavierStokesDisc

	--diffusion and convection
	-- create dirichlet boundary for concentration
	--dirichletBND = DirichletBoundary()
	--dirichletBND:add("ConcentrationDirichletBnd", "c", walls)
	Porosity = 0.1
	rhophi = Porosity

	--Diffusion = ScaleAddLinkerMatrix()
	--Diffusion:add(Porosity, MolecularDiffusion)

	--TransportEq = ConvectionDiffusion("c", "Inner", discType)
	--TransportEq:set_mass_scale(rhophi)
	--TransportEq:set_velocity(NavierStokesDisc)
	--TransportEq:set_diffusion(Diffusion)
	print("Transport Equation created.")
	--
		
	local porder = approxSpace:lfeid(dim):order()
	local vorder = approxSpace:lfeid(0):order()
	local corder = 1 
	
	--upwind if available
	if discType == "fv1" or discType == "fvcr" then
		NavierStokesDisc:set_upwind(upwind)
		NavierStokesDisc:set_peclet_blend(bPecletBlend)
	end
	
	-- fv1 must be stablilized
	if discType == "fv1" then
		NavierStokesDisc:set_stabilization(stab, diffLength)
		NavierStokesDisc:set_pac_upwind(true)
	end
	
	-- fe must be stabilized for (Pk, Pk) space
	if discType == "fe" and porder == vorder then
		NavierStokesDisc:set_stabilization(3)
	end
	if discType == "fe" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	if discType == "fv" then
		NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
	end
	
	-- setup Outlet
	--OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
	--OutletDisc:add("Outlet")
	
	-- setup Inlet
	function inletVel2d(x, y, t)
		return 4 * Um * y * (H-y) / (H*H), 0.0 
	end
	function inletVel3d(x, y, z, t)
		return 16 * Um * y * z * (B-y) * (H-z) / (B*B*H*H), 0.0, 0.0 --return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0 --da quadratisch
	end
	InletDisc = NavierStokesInflow(NavierStokesDisc)
	InletDisc:add("inletVel"..dim.."d", "Inlet, Outlet")
	
	--setup Walles
	WallDisc = NavierStokesWall(NavierStokesDisc)
	WallDisc:add(walls)

	
	-- Finally we create the discretization object which combines all the
	-- separate discretizations into one domain discretization.
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(NavierStokesDisc)
	--domainDisc:add(TransportEq)
	domainDisc:add(InletDisc)
	domainDisc:add(WallDisc)
	--domainDisc:add(OutletDisc)
	
	return domainDisc
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------

function CreateSolver(approxSpace, discType, p)

	local base = LU()
	
	local smoother = nil
	if discType == "fvcr" or discType == "fecr" then 
		smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	elseif discType == "fv1" then 
		smoother = ILU()
		smoother:set_damp(0.7)
	elseif discType == "fe" and porder == vorder then
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
	end
	
	local smooth = util.smooth.parseParams()
	smoother = util.smooth.create(smooth)

	local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
	local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
							 cycle, base, baseLev, bRAP)
	gmg:add_prolongation_post_process(AverageComponent("p"))
	-- transfer = StdTransfer()
	-- transfer:enable_p1_lagrange_optimization(false)
	-- gmg:set_transfer(transfer)

	local sol = "gmres"--"lu" --"bicgstab"-- util.solver.parseParams() --
	local solver = util.solver.create(sol, gmg)
	local eps = 1e-2--4 bei 3d
	if bStokes then
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))--(4, 5*eps, 1e-99, true, true)) --4 steps bei 3d 
	else 
		solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true))--(1, 5*eps, 1e-2, true, true))	
	end
		
	local convCheck = ConvCheck(500, 1e-11, 1e-99, true)--(1, 10*eps, 1e-99, true, true) 
	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)
	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
	AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Run Problem
--------------------------------------------------------------------------------

local p = vorder
local dom = CreateDomain()
local approxSpace = CreateApproxSpace(dom, discType, vorder, porder,corder)
local domainDisc = CreateDomainDisc(approxSpace, discType, p)
local solver = CreateSolver(approxSpace, discType, p)
--solver:set_debug(GridFunctionDebugWriter(approxSpace))


print(solver:config_string())

local u = GridFunction(approxSpace)
u:set(0)

--	ComputeNonLinearSolution(u, CreateDomainDisc(approxSpace, "fe", p), solver)
ComputeNonLinearSolution(u, domainDisc, solver)

local VelCmp, FctCmp
if 		dim == 1 then VelCmp = {"u"}; 			FctCmp = {"u", "p"}
elseif  dim == 2 then VelCmp = {"u", "v"};	 	FctCmp = {"u", "v", "p"}
elseif  dim == 3 then VelCmp = {"u", "v", "w"}; FctCmp = {"u", "v", "w", "p"}
else print("Choosen Dimension " .. dim .. "not supported. Exiting.") exit() end

vtkWriter = VTKOutput()
vtkWriter:select(VelCmp, "velocity")
vtkWriter:select("p", "pressure")
vtkWriter:select("c", "concentration")
vtkWriter:print(gridName, u)

