from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout


def run_water(fname, equilibrationSteps, steps, save_steps, load_checkpoint=None, platform='CPU'):

    temperature = 320*kelvin
    pressure = 1.0*atmospheres

    checkpoint_interval = 5e5

    gro = GromacsGroFile('in/box.gro')
    top = GromacsTopFile('in/topol.top', unitCellDimensions=gro.getUnitCellDimensions(),
                includeDir='/usr/local/gromacs/share/gromacs/top')
            
    # System Configuration
    nonbondedMethod = PME
    nonbondedCutoff = 1.0*nanometers
    ewaldErrorTolerance = 0.0005
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 0.000001
    hydrogenMass = 1.0*amu

    # Integration Options
    dt = 0.002*picoseconds
    friction = 1.0/picosecond

    barostatInterval = 25

    if platform == 'CPU':
        platform = Platform.getPlatformByName('CPU')
        platformProperties = {}
    elif platform == 'CUDA'
        platform = Platform.getPlatformByName('CUDA')
        platformProperties = {'Precision': 'single'}
    else:
        raise NotImplementedError

    # Prepare the Simulation
    print('Building system...')
            
    system = top.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints,
             ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)

    system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

    integrator = NoseHooverIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
            
    simulation = Simulation(top.topology, system, integrator, platform, platformProperties)

    if load_checkpoint is None:
        gro = GromacsGroFile('in/box.gro')
        positions = gro.positions
        simulation.context.setPositions(positions)
        # Minimize and Equilibrate
        print('Performing energy minimization...')
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temperature)
    else:
        simulation.loadCheckpoint(f'out/checkpoint_{load_checkpoint}.chk')


    dcdReporter = DCDReporter('out/trajectory_'+fname+'.dcd', save_steps)
    dataReporter = StateDataReporter('out/log_'+fname+'.txt', save_steps, totalSteps=equilibrationSteps+steps,
                step=True, speed=True, progress=True, totalEnergy=True, potentialEnergy=True, temperature=True, density=True, separator='\t')
    checkpointReporter = CheckpointReporter(f'out/checkpoint_{fname}.chk', checkpoint_interval)

    print('Equilibrating...')
    simulation.reporters.append(checkpointReporter)
    simulation.reporters.append(dataReporter)
    simulation.step(equilibrationSteps)
    # Simulate

    print('Simulating...')
    simulation.reporters.append(dataReporter)

    simulation.reporters.append(dcdReporter)

    simulation.step(steps)


run_water('320K_1', equilibrationSteps=5e5, steps=1e6, save_steps=1000)
# run_water('320K_2', equilibrationSteps=0, steps=2e5, save_steps=100, load_checkpoint='320K_1')
# run_water('320K_3', equilibrationSteps=0, steps=2e4, save_steps=10, load_checkpoint='320K_1')
