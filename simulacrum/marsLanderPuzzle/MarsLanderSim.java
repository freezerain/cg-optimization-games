package simulacrum.marsLanderPuzzle;
import Archive.MarsLander.*;

import marsLander.Lander;
import marsLander.LanderPath;
import src.DrawerController;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.*;

import static java.lang.Thread.sleep;

public class MarsLanderSim {
    private final static File TEST_FILE = new File(
            MarsLanderSim.class.getResource("simulacrum" + "/marsLanderPuzzle/test1.txt")
                    .getPath());
    private static final int TEST_NUMBER = 4;
    private static boolean IS_SIMULATE = true;
    private static boolean IS_DRAW_VISUAL = true;
    private static int DELAY = 100;

    public static void main(String[] args) throws FileNotFoundException, InterruptedException {
        Scanner sc = new Scanner(TEST_FILE);
        for (int i = 0; i < TEST_NUMBER * 2; i++) sc.nextLine();
        double[] landscape = Arrays.stream(sc.nextLine().split(" "))
                .mapToDouble(Double::parseDouble)
                .toArray();

        double[] landingSite = new double[4];
        double lastX = -1;
        double lastY = -1;
        for (int i = 0; i < landscape.length; i += 2){
            if (landscape[i + 1] == lastY && landscape[i] - lastX >= 1000){
                landingSite[0] = lastX+100;
                landingSite[1] = lastY;
                landingSite[2] = landscape[i]-100;
                landingSite[3] = landscape[i + 1];
                break;
            }
            lastX = landscape[i];
            lastY = landscape[i + 1];
        }
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        evaluate(landscape,landingSite, s1);
        //evaluateStat(landscape, landingSite, s1);
        //testWeights(landscape, landingSite, s1);
        System.out.println("Algho completed!");
    }



    private static void testWeights(double[] landscape, double[] landingSite, int[] initialData) throws InterruptedException {
        System.out.println("Start weight tests");
/*        evaluateCustomTest(new double[]{100.0, 5.0, 0.1, 1});
        evaluateCustomTest(new double[]{100.0, 5.0, 10.0, 1});
        evaluateCustomTest(new double[]{100.0, 0.5, 1.0, 1});
        evaluateCustomTest(new double[]{100.0, 50.0, 1.0, 1});
        evaluateCustomTest(new double[]{1000.0, 5.0, 1.0, 1});
        evaluateCustomTest(new double[]{50.0, 5.0, 1.0, 1});
        evaluateCustomTest(new double[]{10.0, 5.0, 1.0, 1});
        evaluateCustomTest(new double[]{1.0, 5.0, 1.0, 1});
        evaluateCustomTest(new double[]{0.1, 5.0, 1.0, 1});*/
        //evaluateCustomTest(new double[]{10.0, 50.0, 0.1, 1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 1.0, 1.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{0.1, 1.0, 1.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{2.0, 1.0, 1.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 0.1, 1.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 2.0, 1.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 1.0, 0.1, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 1.0, 2.0, 1.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 1.0, 1.0, 0.1, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 1.0, 1.0, 2.0, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{0.2, 0.3, 0.4, 0.1, 0.1});
        evaluateCustomTest(landscape, landingSite, initialData, new double[]{1.0, 2.0, 3.0, 0.5, 0.1});
    }

    private static void evaluateCustomTest(double[] landscape, double[] landingSite, int[] initialData,
                                           double[] weights) throws InterruptedException {
        long start = System.currentTimeMillis();
        Archive.MarsLander.Genetic.State state = new Archive.MarsLander.Genetic.State(initialData[0], initialData[1], initialData[2],
                initialData[3], initialData[4], initialData[5], initialData[6], landscape, landingSite);
        int simSize = 10;
        double[] testAverage = new double[6];
        MarsLanderStat test = new MarsLanderStat();
        test.GENE_WEIGHTS = weights;
        ExecutorService es = Executors.newFixedThreadPool(12);
        List<TestCallable> callList = new ArrayList<>();
        for (int j = 0; j < simSize; j++)
             callList.add(new TestCallable(test, state));
        List<Future<double[]>> futures = es.invokeAll(callList);
        es.shutdown();
        for (int j = 0; j < futures.size(); j++){
            try {
                double[] result = futures.get(j).get();
                for (int l = 0; l < 6; l++)
                     testAverage[l] += result[l];
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        for (int j = 0; j < testAverage.length; j++)
             testAverage[j] /= simSize;
        System.out.println("Test finished in " + (System.currentTimeMillis() - start) + ".ms");
        System.out.println("Time: " + testAverage[0] + ", " +
                           "First t: " + testAverage[1] + ", " +
                           "Ten t: " + testAverage[2] + ", " +
                           "Best fit: " + testAverage[3]+
                           ", Last fit"+ testAverage[4]+
                           ", Failed" + testAverage[5]);
    }


    private static void evaluateStat(double[] landscape, double[] landingSite, int[] initialData) throws InterruptedException {
        long start = System.currentTimeMillis();
        int simSize = 10;
        Archive.MarsLander.Genetic.State state = new Archive.MarsLander.Genetic.State(initialData[0], initialData[1], initialData[2],
                initialData[3], initialData[4], initialData[5], initialData[6], landscape, landingSite);
        List<List<MarsLanderStat>> testList = new ArrayList<>();
        //testList.add(MarsLanderStat.TestType.REMOVE_DUPLICATES.getTest(new ArrayList<>()));
        for (MarsLanderStat.TestType type : MarsLanderStat.TestType.values())
            if(type!=MarsLanderStat.TestType.WEIGHTS)testList.add(type.getTest(new ArrayList<>()));
        
        for (int k = 0; k < testList.size(); k++){
            System.out.print("Test: " + k + "/" + (testList.size() - 1) + " -> ");
            List<MarsLanderStat> tests = testList.get(k);
            System.out.print(tests.get(0).type.name() + " size:" + (tests.size()) + " -> ");
            
            for (int i = 0; i < tests.size(); i++){
                System.out.print(i + " ");
                MarsLanderStat test = tests.get(i);
                ExecutorService es = Executors.newFixedThreadPool(12);
                List<TestCallable> callList = new ArrayList<>();
                for (int j = 0; j < simSize; j++)
                     callList.add(new TestCallable(test, state));
                List<Future<double[]>> futures = es.invokeAll(callList);
                es.shutdown();
                double[] sum = new double[6];
                for (int j = 0; j < futures.size(); j++){
                    try {
                        double[] result = futures.get(j).get();
                        for (int l = 0; l < 6; l++){
                            sum[l] += result[l];
                        }
                    } catch (ExecutionException e) {
                        e.printStackTrace();
                    }
                }
                    test.numbers[0] =  sum[0] / futures.size();
                    test.numbers[1] =  sum[1] / futures.size();
                    test.numbers[2] =  sum[2] / futures.size();
                    test.numbers[3] = sum[3] /   futures.size();
                    test.numbers[4] = sum[4] /   futures.size();
                    test.numbers[5] =  sum[5];
            }
            System.out.println();
        }
        System.out.println("Test finished in " + (System.currentTimeMillis() - start) + ".ms");
        printStat(testList);
    }

    private static void printStat(List<List<MarsLanderStat>> tests) {
        DrawerController dc = new DrawerController();

        Map<String, Map<String, Map<Boolean, Number>>> booleanTestMap = new HashMap<>();
        Map<String, Map<String, Map<Number, Number>>> numberTestMap = new HashMap<>();
        Map<String, Map<String, Map<double[], Number>>> arrayTestMap = new HashMap<>();

        for (List<MarsLanderStat> testStat : tests)
            for (MarsLanderStat statData : testStat){
                switch (statData.type) {
                    case RANDOM_CROSSOVER:
                    case REMOVE_DUPLICATES: {
                        statData.putBooleanMap(booleanTestMap);
                        break;
                    }
                    case ELITISM_PERCENTAGE:
                    case DESIRED_POPULATION_SIZE:
                    case INDIVIDIAL_LENGTH:
                    case CROSSOVER_PERCENTAGE:
                    case MUTATION_CHANCE: {
                        statData.putNumberMap(numberTestMap);
                        break;
                    }
                    case WEIGHTS: {
                        statData.putArrayMap(arrayTestMap);
                        break;
                    }
                }
            }
        booleanTestMap.forEach(dc::showGraph);
        numberTestMap.forEach(dc::showLineChart);
        arrayTestMap.forEach(dc::showMultiLineChart);
    }

    private static void evaluate(double[] landscape, double[] landingSite, int[] initialData) throws InterruptedException {
        Archive.MarsLander.Genetic.State state = new Archive.MarsLander.Genetic.State(initialData[0], initialData[1], initialData[2],
                initialData[3], initialData[4], initialData[5], initialData[6], landscape, landingSite);
        DrawerController dc = new DrawerController();
        if (IS_DRAW_VISUAL) dc.startApp(landscape);
        int counter = 0;
        Archive.MarsLander.Genetic.Settings settings = new Archive.MarsLander.Genetic.Settings();
        Archive.MarsLander.Genetic.Solver solver = new Archive.MarsLander.Genetic.Solver(settings);
        int[] crossovers = new int[1];
        sleep(4000);
        Archive.MarsLander.Genetic.Individual solution = null;
        while (!state.isLanded) {
            crossovers[0] = 0;
            List<Archive.MarsLander.Genetic.Individual> population  = solver.evolve(solution, state, System.currentTimeMillis(), crossovers);
            System.out.println("crossovers: " + crossovers[0]);
            Archive.MarsLander.Genetic.Individual bestInd = population.get(0);
            System.out.println("best fit:"+bestInd.fitness);
            System.out.println(bestInd.finalState);
            Archive.MarsLander.Genetic.State startState = new Archive.MarsLander.Genetic.State(state);
            if (IS_SIMULATE) {
                int angle = -state.angle;
                int power = 4;
                if(bestInd.finalState.isSafeLanded || solution!=null) {
                    if(solution==null || bestInd.fitness>solution.fitness)
                        solution = bestInd;
                    angle = solution.genes.get(0).angle;
                    power = solution.genes.get(0).power;
                }
                else{
                    System.out.println("---- default move!!!! ----");
                    if(state.x>landingSite[2]|| state.hS> 50.0) 
                        if(state.hS>-50.0)angle = (int) (state.angle > 21.9? Math.max(-15, 21.9 - state.angle ):
                            Math.max(15, 21.9 - state.angle));
                    else if(state.x<landingSite[0] || state.hS<-50.0)
                            if(state.hS<50.0)
                        angle = (int) (state.angle > -21.9? Math.max(-15, -21.9 - state.angle):
                                Math.max(15, -21.9 - state.angle));
                    else if(state.hS>20.0)
                        angle = (int) (state.angle > 21.9? Math.max(-15, 21.9 - state.angle ):
                                Math.max(15, 21.9 - state.angle));
                    
                    else if(state.hS<-20.0)
                        angle = (int) (state.angle > -21.9? Math.max(-15, -21.9 - state.angle):
                                Math.max(15, -21.9 - state.angle));
                    
                }
                Archive.MarsLander.Genetic.Gene newGene = new Archive.MarsLander.Genetic.Gene(angle, power);
                state.simulate(settings, newGene);
            }
            
            List<LanderPath> parsed = new ArrayList<>();
            for (Archive.MarsLander.Genetic.Individual ind : population){
                ArrayList<Lander> paths = new ArrayList<>();
                Archive.MarsLander.Genetic.State history = new Archive.MarsLander.Genetic.State(startState);
                for (Archive.MarsLander.Genetic.Gene gene : ind.genes){
                    history.simulate(settings, gene);
                    paths.add(new Lander(history.x, history.y, history.hS, history.vS, history.fuel, history.angle, history.power,
                            history.isLanded, history.isSafeLanded));
                }
                parsed.add(new LanderPath(paths, ind.fitness));
            }
            if (solution!=null){
                ArrayList<Lander> paths = new ArrayList<>();
                Archive.MarsLander.Genetic.State history = new Archive.MarsLander.Genetic.State(startState);
                for (Archive.MarsLander.Genetic.Gene gene : solution.genes){
                    history.simulate(settings, gene);
                    paths.add(new Lander(history.x, history.y, history.hS, history.vS, history.fuel, history.angle, history.power,
                            history.isLanded, history.isSafeLanded));
                }
                parsed.add(new LanderPath(paths, solution.fitness));
            }
            if(solution!=null)
                solution.genes.remove(0);
            if (IS_DRAW_VISUAL) dc.updateCoord(counter++, parsed);
            sleep(DELAY);
        }
    }

    static class TestCallable implements Callable<double[]> {
        final MarsLanderStat stat;
        private final Archive.MarsLander.Genetic.State prototype;

        public TestCallable(MarsLanderStat stat, Archive.MarsLander.Genetic.State prototype) {
            this.stat        = stat;
            this.prototype   = prototype;
        }

        @Override
        public double[] call() {
            Archive.MarsLander.Genetic.Settings s = new Archive.MarsLander.Genetic.Settings();
            s.RANDOM_CROSSOVER_ON_DUPLICATE = stat.RANDOM_CROSSOVER_ON_DUPLICATE;
            s.REMOVE_DUPLICATES             = stat.REMOVE_DUPLICATES;
            s.TOURNAMENT_SIZE               = stat.TOURNAMENT_SIZE;
            s.ELITISM_PERCENTAGE            = stat.ELITISM_PERCENTAGE;
            s.DESIRED_POPULATION_SIZE       = stat.DESIRED_POPULATION_SIZE;
            s.INDIVIDUAL_LENGTH             = stat.INDIVIDUAL_LENGTH;
            s.CROSSOVER_PERCENTAGE          = stat.CROSSOVER_PERCENTAGE;
            s.MUTATION_CHANCE               = stat.MUTATION_CHANCE;
            s.GENE_WEIGHTS                  = stat.GENE_WEIGHTS;
            Archive.MarsLander.Genetic.Solver solver = new Archive.MarsLander.Genetic.Solver(s);

            Archive.MarsLander.Genetic.State gs = new Archive.MarsLander.Genetic.State(prototype);
            int counter = 0;
            int firstFoundCounter = -1;
            int tenSolutionCounter = -1;
            double bestFitness = -1;
            double lastFitness = -1;

            int[] crossovers = new int[1];
            Archive.MarsLander.Genetic.Individual solution = null;
            while (counter < 1000) {
                if(solution!=null)
                    solution.genes.remove(0);
                List<Archive.MarsLander.Genetic.Individual> population = solver.evolve(solution, gs,
                        System.currentTimeMillis(), crossovers);
                Archive.MarsLander.Genetic.State bestState = population.get(0).finalState;
                long solutions = population.stream()
                        .filter(ind -> ind.finalState.isSafeLanded)
                        .count();
                if (solutions > 0){
                    lastFitness = population.get(0).fitness;
                    if (lastFitness > bestFitness) bestFitness = lastFitness;
                }
                if (solutions > 0 && firstFoundCounter < 0) firstFoundCounter = crossovers[0];
                if (solutions >= 10){
                    tenSolutionCounter = crossovers[0];
                    break;
                }
                counter++;
            }
            return new double[]{crossovers[0], firstFoundCounter, tenSolutionCounter,
                        bestFitness, lastFitness, tenSolutionCounter==-1?1:0};
        }
    }
}

