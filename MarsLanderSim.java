import marsLander.Lander;
import marsLander.LanderPath;
import simulacrum.marsLanderPuzzle.MarsLanderStat;
import src.DrawerController;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

import static java.lang.Thread.sleep;

public class MarsLanderSim {
    private final static File TEST_FILE = new File(
            MarsLanderSim.class.getResource("simulacrum" + "/marsLanderPuzzle/test1.txt")
                    .getPath());
    private static final int TEST_NUMBER = 4;
    private static boolean IS_SIMULATE = false;
    private static boolean IS_DRAW_VISUAL = false;
    private static int DELAY = 1;

    public static void main(String[] args) throws FileNotFoundException, InterruptedException {
        Scanner sc = new Scanner(TEST_FILE);
        for (int i = 0; i < TEST_NUMBER * 2; i++){
            sc.nextLine();
        }
        double[] landscape = Arrays.stream(sc.nextLine().split(" "))
                .mapToDouble(Double::parseDouble)
                .toArray();
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        Player.GameState gs = new Player.GameState(s1[0], s1[1], s1[2], s1[3], s1[4], s1[5], s1[6]);
        //evaluate(landscape, gs);
        //evaluateStat(landscape, gs);
        printStat(MarsLanderStat.getTests());
        //evaluateAndTest(gs, player);
        System.out.println("Algho completed!");
    }

    private static void evaluateStat(double[] landscape, Player.GameState prototype) throws InterruptedException {
        long start = System.currentTimeMillis();
        List<List<MarsLanderStat>> testList = MarsLanderStat.getTests();
        for (int k = 0; k < testList.size(); k++){
            System.out.print("Testing: " + k + "/" + testList.size() + " -> ");
            List<MarsLanderStat> tests = testList.get(k);
            for (int i = 0; i < tests.size(); i++){
                System.out.print(i + " ");
                MarsLanderStat test = tests.get(i);
                int counterSum = 0;
                int firstSum = 0;
                int tenSum = 0;
                double bestFitSum = 0;
                double lastFitSum = 0;
                int failed = 0;
                for (int j = 0; j < 10; j++){
                    Player player = new Player();
                    player.simulationInit(landscape);
                    player.GE = new Player.GeneticAlgorithm();
                    player.GE.IS_TOURNAMENT_SELECT = test.IS_TOURNAMENT_SELECT;
                    player.GE.IS_POINT_CROSSOVER = test.IS_POINT_CROSSOVER;
                    player.GE.RANDOM_CROSSOVER_ON_DUPLICATE = test.RANDOM_CROSSOVER_ON_DUPLICATE;
                    player.GE.REMOVE_DUPLICATES = test.REMOVE_DUPLICATES;
                    player.GE.TOURNAMENT_SIZE = test.TOURNAMENT_SIZE;
                    player.GE.ELITISM_PERCENTAGE = test.ELITISM_PERCENTAGE / 100.0;
                    player.GE.DESIRED_POPULATION_SIZE = test.DESIRED_POPULATION_SIZE;
                    player.GE.CROSSOVER_PERCENTAGE = test.CROSSOVER_PERCENTAGE / 100.0;
                    player.GE.MUTATION_CHANCE = test.MUTATION_CHANCE / 100.0;
                    player.GE.GENE_WEIGHTS = test.GENE_WEIGHTS;
                    Player.GameState gs = new Player.GameState(prototype);
                    int counter = 0;
                    int[] counterArr = new int[1];
                    int firstFoundCounter = -1;
                    int tenSolutionCounter = -1;
                    double bestFitness = -1;
                    double lastFitness = -1;
                    while (counter < 100) {
                        List<Player.Individual> indList = player.GE.evaluate(gs,
                                System.currentTimeMillis(), counterArr);
                        long solutions = indList.stream()
                                .filter(ind -> ind.gameStateList.get(
                                        ind.gameStateList.size() - 1).isSafeLanded)
                                .count();
                        if (solutions > 0){
                            lastFitness = indList.get(0).fitnessScore;
                            if (lastFitness > bestFitness) bestFitness = lastFitness;
                        }
                        if (solutions > 0 && firstFoundCounter < 0)
                            firstFoundCounter = counterArr[0];
                        if (solutions >= 10){
                            tenSolutionCounter = counterArr[0];
                            break;
                        }
                        counter++;
                    }
                    if (tenSolutionCounter != -1){
                        counterSum += counter;
                        firstSum += firstFoundCounter;
                        tenSum += tenSolutionCounter;
                        bestFitSum += bestFitness;
                        lastFitSum += lastFitness;
                    } else failed++;

                }

                test.time = counterSum / 100 / (10 - failed);
                test.firstSum = firstSum / (10 - failed);
                test.tenSum = tenSum / (10 - failed);
                test.bestFitSum = bestFitSum / (10 - failed);
                test.lastFitSum = lastFitSum / (10 - failed);
                test.failed = failed;
            }
            System.out.println();
        }
        System.out.println("Test finished in " + (System.currentTimeMillis() - start) + ".ms");
        printStat(testList);
    }

    private static void printStat(List<List<MarsLanderStat>> tests) {
        DrawerController dc = new DrawerController();
        dc.startGraphApp();

        
    }

    private static void evaluate(double[] landscape, Player.GameState prototype) throws InterruptedException {
        Player.GameState gs = new Player.GameState(prototype);
        DrawerController dc = new DrawerController();
        Player player = new Player();
        player.simulationInit(landscape);
        dc.startApp(landscape);
        sleep(3000);
        int counter = 0;
        while (!gs.isLanded) {
            List<Player.Individual> indList = player.GE.evaluate(gs, System.currentTimeMillis(),
                    new int[1]);
            if (IS_SIMULATE) gs = indList.get(0).gameStateList.get(0);
            List<LanderPath> parsed = indList.stream()
                    .map(l -> new LanderPath(l.gameStateList.stream()
                            .map(s -> new Lander(s.x, s.y, s.hSpeed, s.vSpeed, s.fuel, s.angle,
                                    s.power, s.isLanded, s.isSafeLanded))
                            .collect(Collectors.toList()), l.fitnessScore))
                    .collect(Collectors.toList());
            if (IS_DRAW_VISUAL) dc.updateCoord(counter++, parsed);
            sleep(DELAY);
        }
    }
}

