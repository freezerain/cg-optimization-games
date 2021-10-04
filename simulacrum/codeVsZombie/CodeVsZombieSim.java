package simulacrum.codeVsZombie;

import simulacrum.marsLanderPuzzle.MarsLanderStat;
import src.DrawerController;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import Archive.CodeVsZombie2;
import Archive.CodeVsZombie2.*;


import static java.lang.Thread.sleep;

public class CodeVsZombieSim {
    //static int[] player = {500, 4500};
/*    static int[] humans = {100,4000,130,5000,10,4500,500,3500,10,5500,100,3000};
    static int[] zombies = {8000,4500,9000,4500,10000,4500,11000,4500,12000,4500,
            13000,4500,14000,4500,15000,3500,14500,2500,15900,500};*/
/*    static int[] player = {0, 4000};
    static int[] humans = {0, 1000, 0, 8000};
    static int[] zombies = {3000, 1000, 3000, 8000, 4000, 1000, 4000, 8000, 5000, 1000, 5000, 8000,
            7000, 1000, 7000, 8000, 9000, 1000, 9000, 8000, 11000, 1000, 11000, 8000, 13000, 1000,
            13000, 8000, 14000, 1000, 14000, 8000, 14500, 1000, 14500, 8000, 15000, 1000, 15000,
            8000};*/

    static int[] player = {3989, 3259};
    static int[] humans = { 647,384,60,1262,1391,1601,1363,422,15470,384,15060,1262,11391,1601,11363,422};
    static int[] zombies = {7900,1579,8500,2470,7500,3798,6500,4682,9000,5664,7500,6319,8500,7094,7800,8447,
            8100,8847,0,7000,1000,7900,3000,8500,5000,7500,7000,6500,9000,7000,11000,7500,13000,8500,15000,7800};


    public static void main(String[] args) throws InterruptedException {
        evaluate();
        // evaluateStat();
        //testWeights();
    }

    private static void evaluate() throws InterruptedException {
        Genetic.Settings settings = new Genetic.Settings();
        Genetic.Solver solver = new Genetic.Solver(settings);
        List<Genetic.Individual> pop = new ArrayList<>();
        int id = 0;
        List<Actor> humansActor = new ArrayList<>();
        for (int i = 0; i < humans.length; i += 2){
            humansActor.add(new Actor(id++, humans[i], humans[i + 1]));
        }
        List<Actor> zombiesActor = new ArrayList<>();
        for (int i = 0; i < zombies.length; i += 2){
            zombiesActor.add(new  Actor(id++, zombies[i], zombies[i + 1]));
        }
        CodeVsZombie2.Actor playerActor = new CodeVsZombie2.Actor(-1, player[0], player[1]);
        DrawerController dc = new DrawerController();
        CodeVsZombie2.Genetic.State startState = new CodeVsZombie2.Genetic.State(playerActor, humansActor,
                zombiesActor);
        while ((!startState.zombies.isEmpty()) && (!startState.humans.isEmpty())) {
            int[] counter = new int[1];

            pop = solver.evolve(pop, startState, System.currentTimeMillis(), counter);
            System.err.println("crossovers: " + counter[0]);
            var scores = new double[pop.size()];
            var paths = new int[pop.size()][];


            for (int i = 0; i < pop.size(); i++){
                CodeVsZombie2.Genetic.Individual individual = pop.get(i);
                scores[i] = individual.state.humans.isEmpty() ? -1.0 : individual.fitness;
                List<CodeVsZombie2.Genetic.Gene> genes = individual.genes;
                int[] path = new int[genes.size() * 2];
                CodeVsZombie2.Genetic.State state = new CodeVsZombie2.Genetic.State(startState);
                for (int j = 0; j < genes.size(); j++){
                    CodeVsZombie2.Genetic.Gene gene = genes.get(j);
                    state.simulate(gene);
                    int nextX = state.player.x;
                    int nextY = state.player.y;
                    path[j * 2]     = nextX;
                    path[j * 2 + 1] = nextY;
                }

                paths[i] = path;
            }


            CodeVsZombie2.Genetic.Individual bestInd = pop.get(0);
            CodeVsZombie2.Genetic.Gene nextGene = pop.get(0).genes.get(0);
            var nextX = Math.max(0, (int) (startState.player.x +
                                           nextGene.distance * Math.cos(nextGene.angle * 22.5)));
            var nextY = Math.max(0, (int) (startState.player.y +
                                           nextGene.distance * Math.sin(nextGene.angle * 22.5)));

            startState = bestInd.startState;
            // System.err.println("stateX/geneX:" + startState.player.x +"/" + nextX+ " 
            // stateY/geneY:" + startState.player.y + "/"+nextY);

            System.err.println("score: " + startState.score);
            humans = new int[startState.humans.size() * 2];
            for (int i = 0; i < startState.humans.size(); i++){
                humans[i * 2]     = startState.humans.get(i).x;
                humans[i * 2 + 1] = startState.humans.get(i).y;
            }
            zombies = new int[startState.zombies.size() * 2];
            for (int i = 0; i < startState.zombies.size(); i++){
                zombies[i * 2]     = startState.zombies.get(i).x;
                zombies[i * 2 + 1] = startState.zombies.get(i).y;
            }

            dc.updateCodeVsZombieApp(new int[]{startState.player.x, startState.player.y}, humans,
                    zombies, paths, scores);
            sleep(1000);
        }
    }

    private static void testWeights() throws InterruptedException {
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
        evaluateCustomTest(new double[]{1.0, 0.5, 0.01, 0.01});
    }

    private static void evaluateCustomTest(double[] weights) throws InterruptedException {
        long start = System.currentTimeMillis();
        int simSize = 100;
        double[] testAverage = new double[4];
        MarsLanderStat test = new MarsLanderStat();
        test.names = new String[]{"humans alive", "crossovers", "game turns", "score"};
        test.GENE_WEIGHTS = weights;
        ExecutorService es = Executors.newFixedThreadPool(12);
        List<TestCallable> callList = new ArrayList<>();
        for (int j = 0; j < simSize; j++)
             callList.add(new TestCallable(test));
        List<Future<double[]>> futures = es.invokeAll(callList);
        es.shutdown();
        for (int j = 0; j < futures.size(); j++){
            try {
                double[] result = futures.get(j).get();
                for (int l = 0; l < 4; l++)
                     testAverage[l] += result[l];
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        for (int j = 0; j < testAverage.length; j++)
             testAverage[j] /= simSize;
        System.out.println("Test finished in " + (System.currentTimeMillis() - start) + ".ms");
        System.out.println("Alive: " + testAverage[0] + ", " +
                           "Crossovers: " + testAverage[1] + ", " +
                           "Turns: " + testAverage[2] + ", " +
                           "Score: " + testAverage[3]);
    }


    private static void evaluateStat() throws InterruptedException {
        long start = System.currentTimeMillis();
        int simSize = 50;

        List<List<MarsLanderStat>> testList = new ArrayList<>();
        testList.add(MarsLanderStat.TestType.REMOVE_DUPLICATES.getTest(new ArrayList<>()));
        testList.add(MarsLanderStat.TestType.RANDOM_CROSSOVER.getTest(new ArrayList<>()));
        // testList.add(MarsLanderStat.TestType.CROSSOVER_PERCENTAGE.getTest(new ArrayList<>()));
        // testList.add(MarsLanderStat.TestType.MUTATION_CHANCE.getTest(new ArrayList<>()));
        //testList.add(MarsLanderStat.TestType.WEIGHTS.getTest(new ArrayList<>()));
        //testList.add(MarsLanderStat.TestType.ELITISM_PERCENTAGE.getTest(new ArrayList<>()));
        //testList.add(MarsLanderStat.TestType.WEIGHTS.getTest(new ArrayList<>()));
        //  for (MarsLanderStat.TestType type : MarsLanderStat.TestType.values())
        //      if(type != MarsLanderStat.TestType.WEIGHTS)testList.add(type.getTest(new 
        //      ArrayList<>()));


        for (int k = 0; k < testList.size(); k++){
            System.out.print("Test: " + k + "/" + (testList.size() - 1) + " -> ");
            List<MarsLanderStat> tests = testList.get(k);
            System.out.print(tests.get(0).type.name() + " size:" + (tests.size()) + " -> ");

            for (int i = 0; i < tests.size(); i++){
                System.out.print(i + " ");
                MarsLanderStat test = tests.get(i);
                test.names = new String[]{"humans alive", "crossovers", "game turns", "score"};
                ExecutorService es = Executors.newFixedThreadPool(12);
                List<CodeVsZombieSim.TestCallable> callList = new ArrayList<>();
                for (int j = 0; j < simSize; j++)
                     callList.add(new CodeVsZombieSim.TestCallable(test));
                List<Future<double[]>> futures = es.invokeAll(callList);
                es.shutdown();
                double[] sum = new double[4];
                for (int j = 0; j < futures.size(); j++){
                    try {
                        double[] result = futures.get(j).get();
                        for (int l = 0; l < 4; l++){
                            sum[l] += result[l];
                        }
                    } catch (ExecutionException e) {
                        e.printStackTrace();
                    }
                }
                test.numbers[0] = sum[0] / futures.size();
                test.numbers[1] = sum[1] / futures.size();
                test.numbers[2] = sum[2] / futures.size();
                test.numbers[3] = sum[3] / futures.size();
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

    static class TestCallable implements Callable<double[]> {
        final MarsLanderStat stat;

        public TestCallable(MarsLanderStat stat) {
            this.stat = stat;
        }

        @Override
        public double[] call() {
            CodeVsZombie2.Genetic.Settings s = new CodeVsZombie2.Genetic.Settings();
            s.RANDOM_CROSSOVER_ON_DUPLICATE = stat.RANDOM_CROSSOVER_ON_DUPLICATE;
            s.REMOVE_DUPLICATES             = stat.REMOVE_DUPLICATES;
            s.TOURNAMENT_SIZE               = stat.TOURNAMENT_SIZE;
            s.ELITISM_PERCENTAGE            = stat.ELITISM_PERCENTAGE;
            s.DESIRED_POPULATION_SIZE       = stat.DESIRED_POPULATION_SIZE;
            s.INDIVIDUAL_LENGTH             = stat.INDIVIDUAL_LENGTH;
            s.CROSSOVER_PERCENTAGE          = stat.CROSSOVER_PERCENTAGE;
            s.MUTATION_CHANCE               = stat.MUTATION_CHANCE;
            s.GENE_WEIGHTS                  = stat.GENE_WEIGHTS;
            var solver = new CodeVsZombie2.Genetic.Solver(s);

            int counter = 0;

            int id = 0;
            List<CodeVsZombie2.Actor> humansActor = new ArrayList<>();
            for (int i = 0; i < humans.length; i += 2)
                 humansActor.add(new CodeVsZombie2.Actor(id++, humans[i], humans[i + 1]));
            List<CodeVsZombie2.Actor> zombiesActor = new ArrayList<>();
            for (int i = 0; i < zombies.length; i += 2)
                 zombiesActor.add(new CodeVsZombie2.Actor(id++, zombies[i], zombies[i + 1]));
            CodeVsZombie2.Actor playerActor = new CodeVsZombie2.Actor(-1, player[0], player[1]);
            CodeVsZombie2.Genetic.State startState = new CodeVsZombie2.Genetic.State(playerActor, humansActor,
                    zombiesActor);
            int[] crossovers = new int[1];
            List<CodeVsZombie2.Genetic.Individual> pop = new ArrayList<>();
            while (startState.humans.size() > 0 && startState.zombies.size() > 0) {
                pop = solver.evolve(pop, startState, System.currentTimeMillis(), crossovers);
                CodeVsZombie2.Genetic.Gene nextGene = pop.get(0).genes.get(0);
                startState.simulate(nextGene);
                counter++;
            }
            return new double[]{startState.humans.size(), crossovers[0], counter, startState.score};
        }
    }
}
