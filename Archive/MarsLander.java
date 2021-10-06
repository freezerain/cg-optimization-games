package Archive;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;

public class MarsLander {

    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        int N = in.nextInt();
        double[] landscape = new double[N * 2];
        double[] landingSite = new double[4];
        boolean isLandFound = false;
        double lastX = -1;
        double lastY = -1;
        for (int i = 0; i < N * 2; i += 2){
            int x = in.nextInt();
            int y = in.nextInt();
            if ((!isLandFound) && y == lastY && x - lastX >= 1000){
                isLandFound    = true;
                landingSite[0] = lastX+50.0;
                landingSite[1] = lastY;
                landingSite[2] = x-50.0;
                landingSite[3] = y;
            }
            lastX            = x;
            lastY            = y;
            landscape[i]     = x;
            landscape[i + 1] = y;
        }
        System.err.println();

        for (int i = 0; i < landscape.length; i++){
            System.err.print(landscape[i] + ", ");
        }
        System.err.println();
        for (int i = 0; i < landingSite.length; i++){
            System.err.print(landingSite[i] + ", ");
        }
        System.err.println();
        
        Genetic.Settings settings = new Genetic.Settings();
        Genetic.Solver solver = new Genetic.Solver(settings);
        Genetic.Individual solution = null;
        int[] crossovers = new int[1];
        while (true) {
            Genetic.State state = new Genetic.State(in.nextInt(), in.nextInt(), in.nextInt(),
                    in.nextInt(), in.nextInt(), in.nextInt(), in.nextInt(), landscape, landingSite);
            List<Genetic.Individual> population = solver.evolve(solution, state, 
                    System.currentTimeMillis(), crossovers);
            Genetic.Individual bestInd = population.get(0);
            System.err.println("pop size: " + population.size()+ "crossovers: "+crossovers[0] +" is win: " +bestInd.finalState.isSafeLanded 
                               + " fitness: " + bestInd.fitness);
            System.err.println(bestInd.finalState);
            System.err.println(solution==null? "no solution" : (solution.finalState + "\nnext gene: " 
                                                                + solution.genes.get(0)));
            crossovers[0]=0;
            
            if(bestInd.finalState.isSafeLanded || solution!=null) {
                int angle = -state.angle;
                int power = 4;
                if(solution==null || bestInd.fitness>solution.fitness)
                    solution = bestInd;
                angle = solution.genes.get(0).angle;
                power = solution.genes.get(0).power;
                state.simulate(settings, new Genetic.Gene(angle, power));
                System.out.println(state.angle + " " + state.power);
            } else {
                System.err.println("---- default move!!!! ----");
                int angle = 0;
                double x = state.x;
                double s = state.hS;
                double lX = landingSite[0];
                double rX = landingSite[2];
                if(x>lX && x<rX)
                    angle = s<-10.0? -15: s>10.0? 15:0;
                else if(Math.abs(s)>30.0)
                    angle = s<-30.0 ? -22:22;
                else angle = x<lX ? -22:22;
                System.out.println(angle +" 4");
            }
            
            if(solution!=null)
                solution.genes.remove(0);
        }
    }
    


    public static class Genetic {
        private static final Random R = new Random();

        public enum CrossoverType {
            CONTINUOUS {
                @Override
                public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                    double crossoverPoint = R.nextDouble();
                    double crossoverRest = 1.0 - crossoverPoint;
                    List<Gene> genes1 = p1.genes;
                    List<Gene> genes2 = p2.genes;
                    List<Gene> newGenes1 = new ArrayList<>(genes1.size());
                    List<Gene> newGenes2 = new ArrayList<>(genes2.size());
                    int lastCrossover = Math.min(genes1.size(), genes2.size());
                    for (int j = 0; j < lastCrossover; j++){
                        Gene gene1 = genes1.get(j);
                        Gene gene2 = genes2.get(j);
                        int angle1 = gene1.angle;
                        int angle2 = gene2.angle;
                        int power1 = gene1.power;
                        int power2 = gene2.power;
                        int newAngle1 = (int) (crossoverPoint * angle1 + crossoverRest * angle2);
                        int newAngle2 = (int) (crossoverPoint * angle2 + crossoverRest * angle1);
                        int newPower1 = (int) (crossoverPoint * power1 + crossoverRest * power2);
                        int newPower2 = (int) (crossoverPoint * power2 + crossoverRest * power1);
                        newGenes1.add(new Gene(newAngle1, newPower1));
                        newGenes2.add(new Gene(newAngle2, newPower2));
                    }
                    childList.add(new Individual(newGenes1));
                    childList.add(new Individual(newGenes2));
                }
            }, POINT {
                @Override
                public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                    List<Gene> genes1 = p1.genes;
                    List<Gene> genes2 = p2.genes;
                    if (genes1.size() < 2 || genes2.size() < 2){
                        childList.add(p1);
                        childList.add(p2);
                        return;
                    }
                    int crossoverIndex = (int) (R.nextDouble() *
                                                Math.min(genes1.size(), genes2.size()));
                    int max = Math.max(genes1.size(), genes2.size());
                    List<Gene> newGenes1 = new ArrayList<>(genes1.subList(0, crossoverIndex));
                    newGenes1.addAll(
                            new ArrayList<>(genes2.subList(crossoverIndex, genes2.size())));
                    List<Gene> newGenes2 = new ArrayList<>(genes2.subList(0, crossoverIndex));
                    newGenes1.addAll(
                            new ArrayList<>(genes1.subList(crossoverIndex, genes1.size())));
                    childList.add(new Individual(newGenes1));
                    childList.add(new Individual(newGenes2));
                }
            };

            public abstract void doCrossover(List<Individual> childList, Individual p1,
                                             Individual p2);
        }

        public enum SelectType {
            TOURNAMENT {
                @Override
                public Individual doSelect(List<Individual> pop, int var) {
                    Individual bestInTournament = null;
                    for (int j = 0; j < var; j++){
                        int index = R.nextInt(pop.size());
                        Individual ind = pop.get(index);
                        if (bestInTournament == null ||
                            ind.fitness > bestInTournament.fitness) bestInTournament = ind;
                    }
                    return bestInTournament;
                }
            }, WHEEL {
                @Override
                public Individual doSelect(List<Individual> pop, int var) {
                    double fitnessSum = 0.0;
                    for (Individual i : pop) fitnessSum += i.fitness;
                    double[] probabilities = new double[pop.size()];
                    double previousProbability = 0.0;
                    for (int i = 0; i < pop.size(); i++){
                        Individual ind = pop.get(i);
                                           previousProbability += ind.fitness / fitnessSum;
                        probabilities[i] = previousProbability;
                    }
                    double selectProbability = Genetic.R.nextDouble();
                    for (int i = 0; i < probabilities.length; i++)
                        if (selectProbability < probabilities[i]) return pop.get(i);
                    return pop.get(0);
                }
            };

            public abstract Individual doSelect(List<Individual> pop, int var);
        }

        public static class State {
            public private final double[] landscape;
            public private final double[] landingSite;
            public double x;
            public double y;
            public double hS;
            public double vS;
            public int fuel;
            public int angle;
            public int power;
            public boolean isLanded = false;
            public boolean isSafeLanded = false;
            public boolean isOutOfBounds = false;

            public State(double x, double y, double hSpeed, double vSpeed, int fuel, int angle,
                         int power, double[] landscape, double[] landingSite) {
                this.x           = x;
                this.y           = y;
                this.hS          = hSpeed;
                this.vS          = vSpeed;
                this.fuel        = fuel;
                this.angle       = angle;
                this.power       = power;
                this.landscape   = landscape;
                this.landingSite = landingSite;
            }

            public State(State gs) {
                this(gs.x, gs.y, gs.hS, gs.vS, gs.fuel, gs.angle, gs.power, gs.landscape,
                        gs.landingSite);
                this.isLanded      = gs.isLanded;
                this.isSafeLanded  = gs.isSafeLanded;
                this.isOutOfBounds = gs.isOutOfBounds;
            }

            public void simulate(Settings s, Gene gene) {
                int dAngle = gene.angle;
                int dPower = gene.power;
                if (dPower != 0) power = Math.max(0, Math.min(4, power + dPower));
                if (dAngle != 0) angle = Math.max(-90, Math.min(90, angle + dAngle));
                double oldX = x;
                double oldY = y;
                move();
                if (isIntersect(oldX, oldY)){
                    isLanded     = true;
                    isSafeLanded = isSafeLand();
                } else if (isOutOfBounds()) isOutOfBounds = true;
            }

            private boolean isOutOfBounds() {
                return x < 0.0 || x > 7000.0 || y < 0.0 || y > 3000.0;
            }

            private boolean isIntersect(double oldX, double oldY) {
                Line2D line = new Line2D.Double(new Point2D.Double(oldX, oldY),
                        new Point2D.Double(x, y));
                for (int i = 0; i < landscape.length - 2; i += 2){
                    if (line.intersectsLine(landscape[i], landscape[i + 1], landscape[i + 2],
                            landscape[i + 3])){
                 /*   double[] intersection = intersection(line,
                            new Line2D.Double(new Point2D.Double(landscape[i], landscape[i + 1]),
                                    new Point2D.Double(landscape[i + 2], landscape[i + 3])));
                    x = intersection[0];
                    y = intersection[1];*/
                        return true;
                    }
                }
                return false;
            }
            
/*
        private double getDistance(double[] inter) {
            double interX = inter[0];
            double interY = inter[1];
            if (interX >= safeAreaBegin && interX <= safeAreaEnd){
                return 0.0;
            }
            double dist = 100.0;
            if (interX < safeAreaBegin){
                for (int i = 0; i < landscape.length; i += 2){
                    if (landscape[i] >= safeAreaBegin) break;
                    if (landscape[i] < interX) continue;
                    if (landscape[i] > interX){
                        dist += Point2D.distance(interX, interY, landscape[i], landscape[i + 1]);
                        interX = landscape[i];
                        interY = landscape[i + 1];
                    }
                }
            } else {
                for (int i = landscape.length - 2; i >= 0; i -= 2){
                    if (landscape[i] <= safeAreaEnd) break;
                    if (landscape[i] > interX) continue;
                    if (landscape[i] < interX){
                        dist += Point2D.distance(interX, interY, landscape[i], landscape[i + 1]);
                        interX = landscape[i];
                        interY = landscape[i + 1];
                    }
                }
            }
            return dist;
        }
*/

            private double[] intersection(Line2D a, Line2D b) {
                double x1 = a.getX1(), y1 = a.getY1(), x2 = a.getX2(), y2 = a.getY2(), x3 =
                        b.getX1(), y3 = b
                        .getY1(), x4 = b.getX2(), y4 = b.getY2();
                double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
                if (d == 0){
                    return new double[]{0.0, 0.0};
                }

                double xi = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
                double yi = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

                return new double[]{xi, yi};
            }

            private boolean isSafeLand() {
                if (x > landingSite[2] || x < landingSite[0] || angle != 0 || Math.abs(vS) > 40 ||
                    Math.abs(hS) > 20 /*|| Math.abs(y-landingSite[3])>0.001*/) return false;
                return true;
            }

            private void move() {
                power = Math.min(power, fuel);
                double radians = Math.toRadians(angle);
                double xAcc = Math.sin(radians) * power;
                double yAcc = (Math.cos(radians) * power - 3.711);
                x += hS - xAcc * 0.5;
                y += vS + yAcc * 0.5;
                hS -= xAcc;
                vS += yAcc;
                fuel -= power;
            }

            @Override
            public String toString() {
                return "State{" + "x=" + x + ", y=" + y + ", hS=" + hS +
                       ", vS=" + vS + ", fuel=" + fuel + ", angle=" + angle + ", power=" + power +
                       ", isLanded=" + isLanded + ", isSafeLanded=" + isSafeLanded +
                       ", isOutOfBounds=" + isOutOfBounds + ", landscape=" + Arrays.toString(landscape) + 
                       ", landingSite=" + Arrays.toString(landingSite) +
                       '}';
            }
        }

        public static class Individual {
            public List<Gene> genes = new ArrayList<>();
            public State startState;
            public State finalState;
            public double fitness;

            public Individual(List<Gene> genes) {
                this.genes.addAll(genes);
            }

            public Individual() {

            }

            public void simulate(State state, Settings settings) {
                startState = state;
                finalState = new State(state);
                int i = 0;
                while (i <= settings.INDIVIDUAL_LENGTH && (!finalState.isLanded) &&
                       (!finalState.isOutOfBounds)) {
                    if (i >= this.genes.size()) genes.add(new Gene());
                    finalState.simulate(settings, genes.get(i));
                    i++;
                }
                if (genes.size() > i) genes = genes.subList(0, i);
                fitness = calculateFitness(settings);
            }

            private double calculateFitness(Settings settings) {
                double distance = 7000.0;
                if (finalState.isLanded) distance = finalState.x < finalState.landingSite[0] ?
                        finalState.landingSite[0] - finalState.x :
                        finalState.x > finalState.landingSite[2] ?
                                finalState.x - finalState.landingSite[2] : 0;
                double hSpeed = Math.abs(finalState.hS);
                double vSpeed = Math.abs(finalState.vS);
                int angle = Math.abs(finalState.angle);
                int fuel = Math.abs(finalState.fuel);

                double normDistance = 10.0 - ((distance / 7000.0) * 10.0);
                //Make score exponential
                normDistance *= normDistance;
                normDistance *= normDistance;
                double normHS = 10.0 - (((hSpeed < 20.0 ? 0 : hSpeed) / 500.0) * 10.0);
                normHS *= normHS;
                normHS *= normHS;
                double normVS = 10.0 - (((vSpeed < 40.0 ? 0 : vSpeed) / 500.0) * 10.0);
                normVS *= normVS;
                normVS *= normVS;
                double normAngle = 10.0 - ((angle / 90.0) * 10.0);
                normAngle *= normAngle;
                normAngle *= normAngle;
                double normFuel = (fuel / (double) startState.fuel) * 10.0;
                normFuel *= normFuel;
                normFuel *= normFuel;
                //Get fitness
                double fitnessScore = 0;
                fitnessScore += normDistance * settings.GENE_WEIGHTS[0];
                fitnessScore += normHS * settings.GENE_WEIGHTS[1];
                fitnessScore += normVS * settings.GENE_WEIGHTS[2];
                fitnessScore += normAngle * settings.GENE_WEIGHTS[3];
                if(finalState.isSafeLanded)  fitnessScore += normFuel * settings.GENE_WEIGHTS[4];
                return fitnessScore;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                Individual that = (Individual) o;

                if (Double.compare(that.fitness, fitness) != 0) return false;
                if (!genes.equals(that.genes)) return false;
                if (startState != null ? !startState.equals(that.startState) :
                        that.startState != null) return false;
                return finalState != null ? finalState.equals(that.finalState) :
                        that.finalState == null;
            }

            @Override
            public int hashCode() {
                int result;
                long temp;
                result = genes.hashCode();
                result = 31 * result + (startState != null ? startState.hashCode() : 0);
                result = 31 * result + (finalState != null ? finalState.hashCode() : 0);
                temp   = Double.doubleToLongBits(fitness);
                result = 31 * result + (int) (temp ^ (temp >>> 32));
                return result;
            }
        }

        public static class Gene {
            public final int angle;
            public final int power;

            public Gene() {
                int angleDir = Genetic.R.nextInt(3);
                //angle = angleDir == 0 ? -15 : angleDir == 1 ? 0 : 15;
                angle = Genetic.R.nextInt(31) - 15;
                power = Genetic.R.nextInt(3) - 1;
            }

            public Gene(int angle, int power) {
                this.angle = angle;
                this.power = power;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                Gene gene = (Gene) o;

                if (angle != gene.angle) return false;
                return power == gene.power;
            }

            @Override
            public int hashCode() {
                int result = angle;
                result = 31 * result + power;
                return result;
            }

            @Override
            public String toString() {
                return "Gene{" + "angle=" + angle + ", power=" + power + '}';
            }
        }

        public static class Solver {
            public Settings settings;

            public Solver(Settings settings) {
                this.settings = settings;
            }

            public List<Individual> evolve(Individual previousBest, State state, long startTime,
                                           int[] crossovers) {
                List<Individual> pop = getRandomIndividuals(settings.DESIRED_POPULATION_SIZE);
                if(previousBest!=null) {
                    List<Gene> previosGenes = new ArrayList<Gene>();
                    for (Gene gene : previousBest.genes)
                        previosGenes.add(new Gene(gene.angle, gene.power));
                    pop.add(new Individual(previosGenes));
                }
                int counter = 0;
                while (System.currentTimeMillis() - startTime < settings.EVALUATE_TIME) {
                    int eliteSize = (int) (pop.size() * settings.ELITISM_PERCENTAGE);
                    List<Individual> elite = new ArrayList<>(pop.subList(0, eliteSize));
                    List<Individual> newPop = crossover(pop);
                    mutate(newPop);
                    elite.addAll(newPop);
                    newPop = elite;
                    newPop = resizeList(newPop);
                    for (Individual i : newPop) i.simulate(state, settings);
                    if (settings.REMOVE_DUPLICATES) newPop = new ArrayList<>(new HashSet<>(newPop));
                    newPop.sort(Comparator.comparingDouble((Individual i) -> i.fitness).reversed());
                    counter++;
                    pop = newPop;
                }
                crossovers[0] += counter;
                return pop;
            }
            

            private void mutate(List<Individual> pop) {
                for (Individual ind : pop){
                    List<Gene> genes = ind.genes;
                    for (int j = 0; j < genes.size(); j++)
                        if (Genetic.R.nextDouble() < settings.MUTATION_CHANCE) genes.set(j,
                                new Gene());
                }
            }
            private List<Individual> crossover(List<Individual> pop) {
                List<Individual> result = new ArrayList<>();
                while (result.size() <
                       settings.DESIRED_POPULATION_SIZE * settings.CROSSOVER_PERCENTAGE) {
                    Individual i1 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                    Individual i2 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                    if (settings.RANDOM_CROSSOVER_ON_DUPLICATE){
                        ArrayList<Gene> randGenes = new ArrayList<>(settings.INDIVIDUAL_LENGTH);
                        for (int i = 0; i < settings.INDIVIDUAL_LENGTH; i++)
                             randGenes.add(new Gene());
                        i2 = new Individual(randGenes);
                    } else {
                        int counter = 0;
                        while (i1 == i2 && counter < 100) {
                            i2 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                            counter++;
                        }
                    }
                    settings.crossoverType.doCrossover(result, i1, i2);
                }
                return result;
            }
            private List<Individual> resizeList(List<Individual> l1) {
                if (l1.size() > settings.DESIRED_POPULATION_SIZE) l1 = new ArrayList<>(
                        l1.subList(0, settings.DESIRED_POPULATION_SIZE));
                else if (l1.size() < settings.DESIRED_POPULATION_SIZE) l1.addAll(
                        getRandomIndividuals(settings.DESIRED_POPULATION_SIZE - l1.size()));
                return l1;
            }

            public List<Individual> getRandomIndividuals(int length) {
                List<Individual> result = new ArrayList<>(length);
                for (int i = 0; i < length; i++)
                     result.add(new Individual());
                return result;
            }
        }

        public static class Settings {
            public SelectType selectType = SelectType.TOURNAMENT;
            public CrossoverType crossoverType = CrossoverType.POINT;
            public boolean RANDOM_CROSSOVER_ON_DUPLICATE = false;
            public boolean REMOVE_DUPLICATES = false;
            public int TOURNAMENT_SIZE = 4;
            public double ELITISM_PERCENTAGE = 0.1;
            public int INDIVIDUAL_LENGTH = 150;
            public int DESIRED_POPULATION_SIZE = 50;
            public double CROSSOVER_PERCENTAGE = 0.7;
            public double MUTATION_CHANCE = 0.01;
            public long EVALUATE_TIME = 95;
            public double[] GENE_WEIGHTS = new double[]{1.0, 2.0, 3.0, 0.5, 0.1};

            public Settings() {
            }

            public Settings(Settings s) {
                selectType                    = s.selectType;
                crossoverType                 = s.crossoverType;
                RANDOM_CROSSOVER_ON_DUPLICATE = s.RANDOM_CROSSOVER_ON_DUPLICATE;
                REMOVE_DUPLICATES             = s.REMOVE_DUPLICATES;
                TOURNAMENT_SIZE               = s.TOURNAMENT_SIZE;
                ELITISM_PERCENTAGE            = s.ELITISM_PERCENTAGE;
                INDIVIDUAL_LENGTH             = s.INDIVIDUAL_LENGTH;
                DESIRED_POPULATION_SIZE       = s.DESIRED_POPULATION_SIZE;
                CROSSOVER_PERCENTAGE          = s.CROSSOVER_PERCENTAGE;
                MUTATION_CHANCE               = s.MUTATION_CHANCE;
                EVALUATE_TIME                 = s.EVALUATE_TIME;
                GENE_WEIGHTS                  = new double[s.GENE_WEIGHTS.length];
                for (int i = 0; i < GENE_WEIGHTS.length; i++)
                     GENE_WEIGHTS[i] = s.GENE_WEIGHTS[i];
            }
        }
    }
}