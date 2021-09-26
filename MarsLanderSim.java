import marsLander.Lander;
import marsLander.LanderPath;
import simulacrum.marsLanderPuzzle.MarsLanderStat;
import src.DrawerController;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

import static java.lang.Thread.sleep;

public class MarsLanderSim {
    private final static File TEST_FILE = new File(
            MarsLanderSim.class.getResource("simulacrum" + "/marsLanderPuzzle/test1.txt")
                    .getPath());

    private static final int TEST_NUMBER = 1;
    private static boolean IS_SIMULATE = false;
    private static boolean IS_DRAW_VISUAL = true;
    private static int DELAY = 500;

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
                landingSite[0] = lastX;
                landingSite[1] = lastY;
                landingSite[2] = landscape[i];
                landingSite[3] = landscape[i + 1];
                break;
            }
            lastX = landscape[i];
            lastY = landscape[i + 1];
        }
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        evaluate(landscape,landingSite, s1);
        //evaluateStat(landscape, landingSite, s1);
        System.out.println("Algho completed!");
    }

    private static void evaluateStat(double[] landscape, double[] landingSite, int[] initialData) throws InterruptedException {
        long start = System.currentTimeMillis();
        int simSize = 3;
        Genetic.State state = new Genetic.State(initialData[0], initialData[1], initialData[2],
                initialData[3], initialData[4], initialData[5], initialData[6]);
        List<List<MarsLanderStat>> testList = new ArrayList<>();
        testList.add(MarsLanderStat.TestType.REMOVE_DUPLICATES.getTest(new ArrayList<>()));
 /*       for (MarsLanderStat.TestType type : MarsLanderStat.TestType.values())
            testList.add(type.getTest(new ArrayList<>()));
        */
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
                     callList.add(new TestCallable(test, landscape, landscape, state));
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
                int failed = (int) sum[5];
                test.numbers[5] = failed;
                if (failed != futures.size()){
                    test.numbers[0] = (int) (sum[0] / (futures.size() - failed));
                    test.numbers[1] = (int) (sum[1] / (futures.size() - failed));
                    test.numbers[2] = (int) (sum[2] / (futures.size() - failed));
                    test.numbers[3] = sum[3] / (futures.size() - failed);
                    test.numbers[4] = sum[4] / (futures.size() - failed);
                }else{
                    test.numbers[0] = 0;
                    test.numbers[1] = 0;
                    test.numbers[2] = 0;
                    test.numbers[3] = 0;
                    test.numbers[4] = 0;
                }
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
        Genetic.State state = new Genetic.State(initialData[0], initialData[1], initialData[2],
                initialData[3], initialData[4], initialData[5], initialData[6]);
        DrawerController dc = new DrawerController();
        if (IS_DRAW_VISUAL) dc.startApp(landscape);
        sleep(3000);
        int counter = 0;
        Genetic.Settings settings = new Genetic.Settings();
        settings.landscape    = landscape;
        settings.landingSite  = landingSite;
        settings.startingFuel = state.fuel;
        Genetic.Solver solver = new Genetic.Solver(settings);
        int[] crossovers = new int[1];
        List<Genetic.Individual> population = solver.getRandomIndividuals(
                settings.DESIRED_POPULATION_SIZE);
        while (!state.isLanded) {
            crossovers[0] = 0;
            population = solver.evolve(population, state, System.currentTimeMillis(), crossovers);
            System.out.println("crossovers: " + crossovers[0]);
            System.out.println(population.get(0).getFitness());
            System.out.println(population.get(0).states.get(population.get(0).states.size() - 1));
            if (IS_SIMULATE) state = population.get(0).states.get(0);
            List<LanderPath> parsed = population.stream()
                    .map(l -> new LanderPath(l.states.stream()
                            .map(s -> new Lander(s.x, s.y, s.hS, s.vS, s.fuel, s.angle, s.power,
                                    s.isLanded, s.isSafeLanded))
                            .collect(Collectors.toList()), l.getFitness()))
                    .collect(Collectors.toList());
            if (IS_DRAW_VISUAL) dc.updateCoord(counter++, parsed);
            sleep(DELAY);
        }
    }

    static class TestCallable implements Callable<double[]> {
        final MarsLanderStat stat;
        private final double[] landscape;
        private final Genetic.State prototype;
        private final double[] landingSite;

        public TestCallable(MarsLanderStat stat, double[] landscape, double[] landingSite,
                            Genetic.State prototype) {
            this.stat        = stat;
            this.landscape   = landscape;
            this.landingSite = landingSite;
            this.prototype   = prototype;
        }

        @Override
        public double[] call() {
            Genetic.Settings s = new Genetic.Settings();
            s.RANDOM_CROSSOVER_ON_DUPLICATE = stat.RANDOM_CROSSOVER_ON_DUPLICATE;
            s.REMOVE_DUPLICATES             = stat.REMOVE_DUPLICATES;
            s.TOURNAMENT_SIZE               = stat.TOURNAMENT_SIZE;
            s.ELITISM_PERCENTAGE            = stat.ELITISM_PERCENTAGE;
            s.DESIRED_POPULATION_SIZE       = stat.DESIRED_POPULATION_SIZE;
            s.CROSSOVER_PERCENTAGE          = stat.CROSSOVER_PERCENTAGE;
            s.MUTATION_CHANCE               = stat.MUTATION_CHANCE;
            s.GENE_WEIGHTS                  = stat.GENE_WEIGHTS;
            s.landscape                     = landscape;
            s.landingSite                   = landingSite;
            s.startingFuel                  = prototype.fuel;
            Genetic.Solver solver = new Genetic.Solver(s);

            Genetic.State gs = new Genetic.State(prototype);
            int counter = 0;
            int firstFoundCounter = -1;
            int tenSolutionCounter = -1;
            double bestFitness = -1;
            double lastFitness = -1;

            int[] crossovers = new int[1];
            List<Genetic.Individual> population = solver.getRandomIndividuals(
                    s.DESIRED_POPULATION_SIZE);
            while (counter < 1000) {
                population = solver.evolve(population, gs,
                        System.currentTimeMillis(), crossovers);
                Genetic.State bestState = population.get(0).states.get(
                        population.get(0).states.size() - 1);
                long solutions = population.stream()
                        .filter(ind -> ind.states.get(ind.states.size() - 1).isSafeLanded)
                        .count();
                if (solutions > 0){
                    lastFitness = population.get(0).getFitness();
                    if (lastFitness > bestFitness) bestFitness = lastFitness;
                }
                if (solutions > 0 && firstFoundCounter < 0) firstFoundCounter = crossovers[0];
                if (solutions >= 10){
                    tenSolutionCounter = crossovers[0];
                    break;
                }
                counter++;
            }
            if (tenSolutionCounter != -1){
                return new double[]{crossovers[0], firstFoundCounter, tenSolutionCounter,
                        bestFitness, lastFitness, 0};
            } else return new double[]{0, 0, 0, 0, 0, 1};
        }
    }
}

