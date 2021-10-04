package simulacrum.marsLanderPuzzle;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MarsLanderStat {
    public TestType type = TestType.DESIRED_POPULATION_SIZE;
    public boolean RANDOM_CROSSOVER_ON_DUPLICATE = true;
    public boolean REMOVE_DUPLICATES = false;
    public int TOURNAMENT_SIZE = 4;
    public double ELITISM_PERCENTAGE = 0.3;
    public int DESIRED_POPULATION_SIZE = 50;
    public int INDIVIDUAL_LENGTH = 20;
    public double CROSSOVER_PERCENTAGE = 0.7;
    public double MUTATION_CHANCE = 0.02;
    public double[] GENE_WEIGHTS = new double[]{100.0, 5.0, 1.0, 1};
    

    public Number[] numbers = new Number[6];
    public String[] names = {"Time", "First time", "Ten time", "Best fitness", "Last fitness", "Failed"};

    public MarsLanderStat() {
    }

    public MarsLanderStat(MarsLanderStat ms) {
        this.type                          = ms.type;
        this.RANDOM_CROSSOVER_ON_DUPLICATE = ms.RANDOM_CROSSOVER_ON_DUPLICATE;
        this.REMOVE_DUPLICATES             = ms.REMOVE_DUPLICATES;
        this.TOURNAMENT_SIZE               = ms.TOURNAMENT_SIZE;
        this.ELITISM_PERCENTAGE            = ms.ELITISM_PERCENTAGE;
        this.DESIRED_POPULATION_SIZE       = ms.DESIRED_POPULATION_SIZE;
        this.CROSSOVER_PERCENTAGE          = ms.CROSSOVER_PERCENTAGE;
        this.MUTATION_CHANCE               = ms.MUTATION_CHANCE;
        this.GENE_WEIGHTS = new double[ms.GENE_WEIGHTS.length];
        for (int i = 0; i < ms.GENE_WEIGHTS.length; i++)
             GENE_WEIGHTS[i] = ms.GENE_WEIGHTS[i];
        this.numbers = new Number[ms.numbers.length];
        for (int i = 0; i < ms.numbers.length; i++)
            numbers[i] = ms.numbers[i];
        this.names = ms.names;
    }

    public MarsLanderStat(TestType type, boolean RANDOM_CROSSOVER_ON_DUPLICATE, boolean REMOVE_DUPLICATES, int TOURNAMENT_SIZE, double ELITISM_PERCENTAGE, int DESIRED_POPULATION_SIZE, double CROSSOVER_PERCENTAGE, double MUTATION_CHANCE, double[] GENE_WEIGHTS) {
        this.type = type;
        this.RANDOM_CROSSOVER_ON_DUPLICATE = RANDOM_CROSSOVER_ON_DUPLICATE;
        this.REMOVE_DUPLICATES = REMOVE_DUPLICATES;
        this.TOURNAMENT_SIZE = TOURNAMENT_SIZE;
        this.ELITISM_PERCENTAGE = ELITISM_PERCENTAGE;
        this.DESIRED_POPULATION_SIZE = DESIRED_POPULATION_SIZE;
        this.CROSSOVER_PERCENTAGE = CROSSOVER_PERCENTAGE;
        this.MUTATION_CHANCE = MUTATION_CHANCE;
        this.GENE_WEIGHTS = GENE_WEIGHTS;
    }
    public static List<MarsLanderStat> getRandomDuplTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            MarsLanderStat newStatTrue = new MarsLanderStat(initStat);
            MarsLanderStat newStatFalse= new MarsLanderStat(initStat);
            newStatTrue.type = TestType.RANDOM_CROSSOVER;
            newStatFalse.type = TestType.RANDOM_CROSSOVER;
            newStatTrue.RANDOM_CROSSOVER_ON_DUPLICATE = true;
            newStatFalse.RANDOM_CROSSOVER_ON_DUPLICATE = false;
            result.add(newStatTrue);
            result.add(newStatFalse);
        }
        return result;
    }
    public static List<MarsLanderStat> getRemoveDuplTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            MarsLanderStat newStatTrue = new MarsLanderStat(initStat);
            MarsLanderStat newStatFalse= new MarsLanderStat(initStat);
            newStatTrue.type = TestType.REMOVE_DUPLICATES;
            newStatFalse.type = TestType.REMOVE_DUPLICATES;
            newStatTrue.REMOVE_DUPLICATES = true;
            newStatFalse.REMOVE_DUPLICATES = false;
            result.add(newStatTrue);
            result.add(newStatFalse);
        }
        return result;
    }
    public static List<MarsLanderStat> getElitismChanceTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            for (int i = 0; i <= 35 ; i+=5){
                MarsLanderStat newStat = new MarsLanderStat(initStat);
                newStat.type = TestType.ELITISM_PERCENTAGE;
                newStat.ELITISM_PERCENTAGE = i/100.0;
                result.add(newStat);
            }
        }
        return result;
    }
    public static List<MarsLanderStat> getPopSizeTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            for (int i = 25; i <= 250 ; i+=25){
                MarsLanderStat newStat = new MarsLanderStat(initStat);
                newStat.type = TestType.DESIRED_POPULATION_SIZE;
                newStat.DESIRED_POPULATION_SIZE = i;
                result.add(newStat);
            }
        }
        return result;
    }
    
    public static List<MarsLanderStat> getIndLengthTest(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            for (int i = 10; i <= 200 ; i+=10){
                MarsLanderStat newStat = new MarsLanderStat(initStat);
                newStat.type = TestType.INDIVIDIAL_LENGTH;
                newStat.INDIVIDUAL_LENGTH = i;
                result.add(newStat);
            }
        }
        return result;
    }
    
    public static List<MarsLanderStat> getCrossoverChanceTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            for (int i = 50; i <= 100-initStat.ELITISM_PERCENTAGE*100 ; i+=5){
                MarsLanderStat newStat = new MarsLanderStat(initStat);
                newStat.type = TestType.CROSSOVER_PERCENTAGE;
                newStat.CROSSOVER_PERCENTAGE = i/100.0;
                result.add(newStat);
            }
        }
        return result;
    }
    public static List<MarsLanderStat> getMutationChanceTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat : startingList){
            for (int i = 1; i <= 5 ; i++){
                MarsLanderStat newStat = new MarsLanderStat(initStat);
                newStat.type = TestType.MUTATION_CHANCE;
                newStat.MUTATION_CHANCE = i/100.0;
                result.add(newStat);
            }
        }
        return result;
    }
    public static List<MarsLanderStat> getWeightTests(List<MarsLanderStat> startingList){
        if(startingList.isEmpty()) startingList.add(new MarsLanderStat());
        ArrayList<MarsLanderStat> result = new ArrayList<>();
        for (MarsLanderStat initStat: startingList) {
            for (int i = 10; i <= 100; i+=10){
                for (int j = 0; j < 5; j++){
                    MarsLanderStat newStat = new MarsLanderStat(initStat);
                    newStat.type = TestType.WEIGHTS;
                    newStat.GENE_WEIGHTS[j] = i/100.0;
                    result.add(newStat);
                }
            }
        }
        return result;
    }
    
    
    
    public void putBooleanMap(Map<String, Map<String, Map<Boolean, Number>>> booleanTestMap){
        Map<String, Map<Boolean, Number>> varMap = booleanTestMap.getOrDefault(
                type.name(), new HashMap<>());
        boolean state = (boolean) type.getTestVar(this);
        for (int i = 0; i < names.length; i++){
            Map<Boolean, Number> testMap = varMap.getOrDefault(names[i],
                    new HashMap<>());
            testMap.put(state, numbers[i]);
            varMap.put(names[i], testMap);
        }
        booleanTestMap.put(type.name(), varMap);
    }
    public void putNumberMap(Map<String, Map<String, Map<Number, Number>>> numberTestMap){
        Map<String, Map<Number, Number>> varMap = numberTestMap.getOrDefault(
                type.name(), new HashMap<>());
        Number n = (Number) type.getTestVar(this);
        for (int i = 0; i < names.length; i++){
            Map<Number, Number> testMap = varMap.getOrDefault(names[i],
                    new HashMap<>());
            testMap.put(n, numbers[i]);
            varMap.put(names[i], testMap);
        }
        numberTestMap.put(type.name(), varMap);
    }
    public void putArrayMap(Map<String, Map<String, Map<double[], Number>>> numberTestMap){
        Map<String, Map<double[], Number>> varMap = numberTestMap.getOrDefault(
                type.name(), new HashMap<>());
        double[] arr = (double[]) type.getTestVar(this);
        for (int i = 0; i < names.length; i++){
            Map<double[], Number> testMap = varMap.getOrDefault(names[i], new HashMap<>());
                testMap.put(arr, numbers[i]);
            varMap.put(names[i], testMap);
        }
        numberTestMap.put(type.name(), varMap);
    }
    
    public enum TestType{
        RANDOM_CROSSOVER{
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getRandomDuplTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.RANDOM_CROSSOVER_ON_DUPLICATE;
            }

        }, 
        REMOVE_DUPLICATES {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getRemoveDuplTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.REMOVE_DUPLICATES;
            }

        }, 
        ELITISM_PERCENTAGE {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getElitismChanceTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.ELITISM_PERCENTAGE;
            }

        }, DESIRED_POPULATION_SIZE {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getPopSizeTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.DESIRED_POPULATION_SIZE;
            }


        }, INDIVIDIAL_LENGTH{
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getIndLengthTest(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.INDIVIDUAL_LENGTH;
            }
        },CROSSOVER_PERCENTAGE {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getCrossoverChanceTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.CROSSOVER_PERCENTAGE;
            }
        }, MUTATION_CHANCE {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getMutationChanceTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.MUTATION_CHANCE;
            }


        },
        WEIGHTS {
            @Override
            public List<MarsLanderStat> getTest(List<MarsLanderStat> startingList) {
                return getWeightTests(startingList);
            }

            @Override
            public Object getTestVar(MarsLanderStat ms) {
                return ms.GENE_WEIGHTS;
            }


        };
        
        
        public abstract List<MarsLanderStat> getTest(List<MarsLanderStat> startingList);
        public abstract Object getTestVar(MarsLanderStat ms);
    }
}
