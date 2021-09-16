package simulacrum.marsLanderPuzzle;

import java.util.ArrayList;
import java.util.List;

public class MarsLanderStat {
    public boolean IS_TOURNAMENT_SELECT = false;
    public boolean IS_POINT_CROSSOVER = false;
    public boolean RANDOM_CROSSOVER_ON_DUPLICATE = true;
    public boolean REMOVE_DUPLICATES = true;
    public int TOURNAMENT_SIZE = 5;
    public double ELITISM_PERCENTAGE = 0.2;
    public int DESIRED_POPULATION_SIZE = 100;
    public double CROSSOVER_PERCENTAGE = 0.6;
    public double MUTATION_CHANCE = 0.02;
    public double[] GENE_WEIGHTS = new double[]{0.3, 0.25, 0.25, 0.2, 1};

    public int time = 0;
    public int firstSum = 0;
    public int tenSum = 0;
    public double bestFitSum = 0;
    public double lastFitSum = 0;
    public int failed = 0;

    public MarsLanderStat() {
    }

    public MarsLanderStat(boolean IS_TOURNAMENT_SELECT, boolean IS_POINT_CROSSOVER,
                          boolean RANDOM_CROSSOVER_ON_DUPLICATE, boolean REMOVE_DUPLICATES,
                          int TOURNAMENT_SIZE, double ELITISM_PERCENTAGE,
                          int DESIRED_POPULATION_SIZE, double CROSSOVER_PERCENTAGE,
                          double MUTATION_CHANCE, double[] GENE_WEIGHTS) {
        this.IS_TOURNAMENT_SELECT = IS_TOURNAMENT_SELECT;
        this.IS_POINT_CROSSOVER = IS_POINT_CROSSOVER;
        this.RANDOM_CROSSOVER_ON_DUPLICATE = RANDOM_CROSSOVER_ON_DUPLICATE;
        this.REMOVE_DUPLICATES = REMOVE_DUPLICATES;
        this.TOURNAMENT_SIZE = TOURNAMENT_SIZE;
        this.ELITISM_PERCENTAGE = ELITISM_PERCENTAGE;
        this.DESIRED_POPULATION_SIZE = DESIRED_POPULATION_SIZE;
        this.CROSSOVER_PERCENTAGE = CROSSOVER_PERCENTAGE;
        this.MUTATION_CHANCE = MUTATION_CHANCE;
        this.GENE_WEIGHTS = GENE_WEIGHTS;
    }

    public static List<List<MarsLanderStat>> getTests() {
        List<List<MarsLanderStat>> testList = new ArrayList<>();
        List<MarsLanderStat> isTourList = new ArrayList<>();
        for (int i = 0; i < 2; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.IS_TOURNAMENT_SELECT = i == 0;
            isTourList.add(m);
        }
        testList.add(isTourList);
        List<MarsLanderStat> isPointList = new ArrayList<>();
        for (int i = 0; i < 2; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.IS_POINT_CROSSOVER = i == 0;
            isPointList.add(m);
        }
        testList.add(isPointList);
        List<MarsLanderStat> isRandList = new ArrayList<>();
        for (int i = 0; i < 2; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.RANDOM_CROSSOVER_ON_DUPLICATE = i == 0;
            isRandList.add(m);
        }
        testList.add(isRandList);
        List<MarsLanderStat> isRemoveList = new ArrayList<>();
        for (int i = 0; i < 2; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.REMOVE_DUPLICATES = i == 0;
            isRemoveList.add(m);
        }
        testList.add(isRemoveList);
        List<MarsLanderStat> tournamentList = new ArrayList<>();
        for (int i = 2; i <= 10; i += 2){
            MarsLanderStat m = new MarsLanderStat();
            m.IS_TOURNAMENT_SELECT = true;
            m.TOURNAMENT_SIZE = i;
            tournamentList.add(m);
        }
        testList.add(tournamentList);
        List<MarsLanderStat> eliteList = new ArrayList<>();
        for (int i = 1; i <= 5; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.ELITISM_PERCENTAGE = i / 100.0;
            eliteList.add(m);
        }
        testList.add(eliteList);
        List<MarsLanderStat> desiredList = new ArrayList<>();
        for (int i = 50; i <= 500; i += 50){
            MarsLanderStat m = new MarsLanderStat();
            m.DESIRED_POPULATION_SIZE = i;
            desiredList.add(m);
        }
        testList.add(desiredList);
        List<MarsLanderStat> crossoverList = new ArrayList<>();
        for (int i = 50; i <= 80; i += 10){
            MarsLanderStat m = new MarsLanderStat();
            m.CROSSOVER_PERCENTAGE = i / 100.0;
            crossoverList.add(m);
        }
        testList.add(crossoverList);
        List<MarsLanderStat> mutationList = new ArrayList<>();
        for (int i = 1; i <= 5; i++){
            MarsLanderStat m = new MarsLanderStat();
            m.MUTATION_CHANCE = i / 100.0;
            mutationList.add(m);
        }
        testList.add(mutationList);
        List<MarsLanderStat> weightList = new ArrayList<>();
        for (int w1 = 10; w1 <= 60; w1 += 10)
            for (int w2 = 10; w2 <= 70 - w1; w2 += 10)
                for (int w3 = 10; w3 <= 80 - w1 - w2; w3 += 10)
                    for (int w4 = 10; w4 <= 90 - w1 - w2 - w3; w4 += 10)
                        for (int w5 = 10; w5 <= 100 - w1 - w2 - w3 - w4; w5 += 10){
                            MarsLanderStat m = new MarsLanderStat();
                            m.GENE_WEIGHTS = new double[]{w1 / 100.0, w2 / 100.0, w3 / 100.0,
                                    w4 / 100.0, w5 / 100.0,};
                            weightList.add(m);
                        }
        testList.add(weightList);
        return testList;
    }

}
