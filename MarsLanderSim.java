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

    private static final int TEST_NUMBER = 0;
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
        //evaluateStat(landscape, gs);
        //printStat(MarsLanderStat.getTests());
        //evaluateAndTest(gs, player);
        System.out.println("Algho completed!");
    }


    private static void evaluateStat(double[] landscape, Player.GameState prototype) throws InterruptedException {
        long start = System.currentTimeMillis();
        int simSize = 10;
        List<List<MarsLanderStat>> testList = MarsLanderStat.getTests();
        for (int k = 0; k < testList.size(); k++){
            System.out.print("Testing: " + k + "/" + (testList.size() - 1) + " -> ");
            List<MarsLanderStat> tests = testList.get(k);
            System.out.print("size: " + (tests.size() - 1));
            for (int i = 0; i < tests.size(); i++){
                System.out.print(i + " ");
                MarsLanderStat test = tests.get(i);
                ExecutorService es = Executors.newFixedThreadPool(12);
                List<TestCallable> callList = new ArrayList<>();
                for (int j = 0; j < simSize; j++){
                    callList.add(new TestCallable(test, landscape, prototype));
                }
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
                test.failed = (int) sum[5];
                if (test.failed != futures.size()){
                    test.time       = (int) (sum[0] / 100 / (futures.size() - test.failed));
                    test.firstSum   = (int) (sum[1] / (futures.size() - test.failed));
                    test.tenSum     = (int) (sum[2] / (futures.size() - test.failed));
                    test.bestFitSum = sum[3] / (futures.size() - test.failed);
                    test.lastFitSum = sum[4] / (futures.size() - test.failed);
                }
            }
            System.out.println();
        }
        System.out.println("Test finished in " + (System.currentTimeMillis() - start) + ".ms");
        printStat(testList);
    }

    private static void printStat(List<List<MarsLanderStat>> tests) {
        DrawerController dc = new DrawerController();
        List<MarsLanderStat> ms = tests.get(0);

        Map<String, Map<String, Map<Boolean, Number>>> booleanTestMap = new HashMap<>();
        Map<String, Map<String, Map<Number, Number>>> numberTestMap = new HashMap<>();
        Map<String, Map<String, Map<double[], Number>>> arrayTestMap = new HashMap<>();

        for (List<MarsLanderStat> testStat : tests)
            for (MarsLanderStat statData : testStat){
                Number[] numbers = new Number[]{statData.time, statData.firstSum, statData.tenSum,
                        statData.bestFitSum, statData.lastFitSum, statData.failed};
                String[] names = {"Time", "First time", "Ten time", "Best fitness", "Last fitness",
                        "Failed"};
                switch (statData.type) {
                    case TOURNAMENT: {
                        boolean state = statData.IS_TOURNAMENT_SELECT;

                        Map<String, Map<Boolean, Number>> varMap = booleanTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Boolean, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(state, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        booleanTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case POINT_CROSSOVER: {
                        boolean state = statData.IS_POINT_CROSSOVER;

                        Map<String, Map<Boolean, Number>> varMap = booleanTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Boolean, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(state, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        booleanTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case RANDOM_CROSSOVER: {
                        boolean state = statData.RANDOM_CROSSOVER_ON_DUPLICATE;

                        Map<String, Map<Boolean, Number>> varMap = booleanTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Boolean, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(state, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        booleanTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case REMOVE_DUPLICATES: {
                        boolean state = statData.REMOVE_DUPLICATES;

                        Map<String, Map<Boolean, Number>> varMap = booleanTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Boolean, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(state, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        booleanTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case TOURNAMENT_SIZE: {
                        Number n = statData.TOURNAMENT_SIZE;

                        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(n, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        numberTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case ELITISM_PERCENTAGE: {
                        Number n = statData.ELITISM_PERCENTAGE;

                        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(n, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        numberTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case DESIRED_POPULATION_SIZE: {
                        Number n = statData.DESIRED_POPULATION_SIZE;

                        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(n, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        numberTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case CROSSOVER_PERCENTAGE: {
                        Number n = statData.CROSSOVER_PERCENTAGE;

                        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(n, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        numberTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    case MUTATION_CHANCE: {
                        Number n = statData.MUTATION_CHANCE;

                        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(n, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        numberTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                    //Fix Weights
                    case WEIGHTS: {
                        double[] arr = statData.GENE_WEIGHTS;

                        Map<String, Map<double[], Number>> varMap = arrayTestMap.getOrDefault(
                                statData.type.name(), new HashMap<>());
                        for (int i = 0; i < names.length; i++){
                            Map<double[], Number> testMap = varMap.getOrDefault(names[i],
                                    new HashMap<>());
                            testMap.put(arr, numbers[i]);
                            varMap.put(names[i], testMap);
                        }
                        arrayTestMap.put(statData.type.name(), varMap);
                        break;
                    }
                }
            }


        //    booleanTestMap.forEach(dc::showGraph);
        numberTestMap.forEach(dc::showLineChart);
        arrayTestMap.forEach(dc::showMultiLineChart);

    }

    private static void evaluate(double[] landscape, double[] landingSite, int[] initialData) throws InterruptedException {
        Genetic.State state = new Genetic.State(initialData[0], initialData[1],
                initialData[2], initialData[3], initialData[4], initialData[5], initialData[6]);
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
        List<Genetic.Individual> population = solver.getRandomIndividuals(settings.DESIRED_POPULATION_SIZE);
        while (!state.isLanded) {
            
            
             population = solver.evolve(population, state,
                    System.currentTimeMillis(), crossovers);
            System.out.println("crossovers: " + crossovers[0]);
            System.out.println(population.get(0));
            if (IS_SIMULATE) state = population.get(0).states.get(0);
            List<LanderPath> parsed = population.stream()
                    .map(l -> new LanderPath(l.states.stream()
                            .map(s -> new Lander(s.x, s.y, s.hS, s.vS, s.fuel, s.angle,
                                    s.power, s.isLanded, s.isSafeLanded))
                            .collect(Collectors.toList()), l.getFitness()))
                    .collect(Collectors.toList());
            if (IS_DRAW_VISUAL) dc.updateCoord(counter++, parsed);
            sleep(DELAY);
        }
    }

    static class TestCallable implements Callable<double[]> {
        final MarsLanderStat stat;
        private double[] landscape;
        private Player.GameState prototype;

        public TestCallable(MarsLanderStat stat, double[] landscape, Player.GameState prototype) {
            this.stat      = stat;
            this.landscape = landscape;
            this.prototype = prototype;
        }

        @Override
        public double[] call() {
            Player player = new Player();
            player.simulationInit(landscape);
            player.GE.IS_TOURNAMENT_SELECT          = stat.IS_TOURNAMENT_SELECT;
            player.GE.IS_POINT_CROSSOVER            = stat.IS_POINT_CROSSOVER;
            player.GE.RANDOM_CROSSOVER_ON_DUPLICATE = stat.RANDOM_CROSSOVER_ON_DUPLICATE;
            player.GE.REMOVE_DUPLICATES             = stat.REMOVE_DUPLICATES;
            player.GE.TOURNAMENT_SIZE               = stat.TOURNAMENT_SIZE;
            player.GE.ELITISM_PERCENTAGE            = stat.ELITISM_PERCENTAGE;
            player.GE.DESIRED_POPULATION_SIZE       = stat.DESIRED_POPULATION_SIZE;
            player.GE.CROSSOVER_PERCENTAGE          = stat.CROSSOVER_PERCENTAGE;
            player.GE.MUTATION_CHANCE               = stat.MUTATION_CHANCE;
            player.GE.GENE_WEIGHTS                  = stat.GENE_WEIGHTS;
            Player.GameState gs = new Player.GameState(prototype);
            int counter = 0;
            int[] counterArr = new int[1];
            int firstFoundCounter = -1;
            int tenSolutionCounter = -1;
            double bestFitness = -1;
            double lastFitness = -1;
            while (counter < 100) {
                List<Player.Individual> indList = player.GE.evaluate(gs, System.currentTimeMillis(),
                        counterArr);
                Player.GameState bestState = indList.get(0).gameStateList.get(
                        indList.get(0).gameStateList.size() - 1);
                long solutions = indList.stream()
                        .filter(ind -> ind.gameStateList.get(
                                ind.gameStateList.size() - 1).isSafeLanded)
                        .count();
                if (solutions > 0){
                    lastFitness = indList.get(0).fitnessScore;
                    if (lastFitness > bestFitness) bestFitness = lastFitness;
                }
                if (solutions > 0 && firstFoundCounter < 0) firstFoundCounter = counterArr[0];
                if (solutions >= 10){
                    tenSolutionCounter = counterArr[0];
                    break;
                }
                counter++;
            }
            if (tenSolutionCounter != -1){
                return new double[]{counterArr[0], firstFoundCounter, tenSolutionCounter,
                        bestFitness, lastFitness, 0};
            } else return new double[]{0, 0, 0, 0, 0, 1};
        }

    }
}

