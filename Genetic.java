import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.*;

public class Genetic {
    private static final Random R = new Random();

    public enum CrossoverType {
        CONTINUOUS {
            @Override
            public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                double crossoverPoint = R.nextDouble();
                double crossoverRest = 1.0 - crossoverPoint;
                List<Gene> genes1 = p1.getGenes();
                List<Gene> genes2 = p2.getGenes();
                List<Gene> newGenes1 = new ArrayList<>(genes1.size());
                List<Gene> newGenes2 = new ArrayList<>(genes2.size());
                int lastCrossover = Math.min(genes1.size(), genes2.size());
                for (int j = 0; j < lastCrossover; j++){
                    Gene gene1 = genes1.get(j);
                    Gene gene2 = genes2.get(j);
                    int[] vars1 = gene1.getValues();
                    int[] vars2 = gene2.getValues();
                    int[] newVars1 = new int[vars1.length];
                    int[] newVars2 = new int[vars1.length];
                    for (int i = 0; i < vars1.length; i++){
                        newVars1[i] = (int) (crossoverPoint * vars1[i] + crossoverRest * vars2[i]);
                        newVars2[i] = (int) (crossoverRest * vars1[i] + crossoverPoint * vars2[i]);
                    }
                    newGenes1.add(new Gene(newVars1));
                    newGenes2.add(new Gene(newVars2));
                }
                childList.add(new Individual(newGenes1));
                childList.add(new Individual(newGenes2));
            }
        }, POINT {
            @Override
            public void doCrossover(List<Individual> childList, Individual p1, Individual p2) {
                List<Gene> genes1 = p1.getGenes();
                List<Gene> genes2 = p2.getGenes();
                if (genes1.size() < 2 || genes2.size() < 2){
                    childList.add(p1);
                    childList.add(p2);
                    return;
                }
                int crossoverIndex = (int) (R.nextDouble() * Math.min(genes1.size(), genes2.size()));
                int max = Math.max(genes1.size(), genes2.size());
                List<Gene> newGenes1 = new ArrayList<>(genes1.subList(0, crossoverIndex));
                newGenes1.addAll(new ArrayList<>(genes2.subList(crossoverIndex, genes2.size())));
                List<Gene> newGenes2 = new ArrayList<>(genes2.subList(0, crossoverIndex));
                newGenes1.addAll(new ArrayList<>(genes1.subList(crossoverIndex, genes1.size())));
                childList.add(new Individual(newGenes1));
                childList.add(new Individual(newGenes2));
            }
        };

        public abstract void doCrossover(List<Individual> childList, Individual p1, Individual p2);
    }

    public static enum SelectType {
        TOURNAMENT {
            @Override
            public Individual doSelect(List<Individual> pop, int var) {
                Individual bestInTournament = null;
                for (int j = 0; j < var; j++){
                    int index = R.nextInt(pop.size());
                    Individual ind = pop.get(index);
                    if (bestInTournament == null ||
                        ind.getFitness() > bestInTournament.getFitness()) bestInTournament = ind;
                }
                return bestInTournament;
            }
        }, WHEEL {
            @Override
            public Individual doSelect(List<Individual> pop, int var) {
                double fitnessSum = 0.0;
                for (Individual i : pop) fitnessSum += i.getFitness();
                double[] probabilities = new double[pop.size()];
                double previousProbability = 0.0;
                for (int i = 0; i < pop.size(); i++){
                    Individual ind = pop.get(i);
                                       previousProbability += ind.getFitness() / fitnessSum;
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
        public double distanceToLanding = Double.MAX_VALUE;
        double x;
        double y;
        double hS;
        double vS;
        int fuel;
        int angle;
        int power;
        boolean isLanded = false;
        boolean isSafeLanded = false;
        boolean isOutOfBounds = false;

        public State(double x, double y, double hSpeed, double vSpeed, int fuel, int angle,
                         int power) {
            this.x             = x;
            this.y             = y;
            this.hS        = hSpeed;
            this.vS        = vSpeed;
            this.fuel          = fuel;
            this.angle         = angle;
            this.power         = power;
        }

        public State(State gs) {
            x                      = gs.x;
            y                      = gs.y;
            hS                = gs.hS;
            vS               = gs.vS;
            fuel                   = gs.fuel;
            angle                  = gs.angle;
            power                  = gs.power;
            isLanded               = gs.isLanded;
            isSafeLanded           = gs.isSafeLanded;
            distanceToLanding = gs.distanceToLanding;
            isOutOfBounds = gs.isOutOfBounds;
        }

        public State simulate(Settings s, Gene gene) {
            State newState = new State(this);
            int dAngle = gene.getValues()[0];
            int dPower = gene.getValues()[1];
            if (dPower != 0) newState.power = Math.max(0, Math.min(4, power + dPower));
            if (dAngle != 0) newState.angle = Math.max(-90, Math.min(90, angle + dAngle));

            int counter = gene.getValues()[2];
            while(counter>0 && (!newState.isLanded) && (!newState.isOutOfBounds)) {
                counter--;
                double oldX = newState.x;
                double oldY = newState.y;
                newState.move();
                if (newState.isIntersect(s, oldX, oldY)){
                    newState.isLanded     = true;
                    newState.isSafeLanded = newState.isSafeLand(s);
                }else if(newState.isOutOfBounds()) newState.isOutOfBounds = true;
            }
            if(counter!=0) gene.getValues()[2] -=counter;
            return newState;
        }

        private boolean isOutOfBounds() {
            return x < 0.0 || x > 7000.0 || y < 0.0 || y > 3000.0;
        }

        private boolean isIntersect(Settings settings, double oldX, double oldY) {
            Line2D line = new Line2D.Double(new Point2D.Double(oldX, oldY),
                    new Point2D.Double(x, y));
            for (int i = 0; i < settings.landscape.length - 2; i += 2){
                if (line.intersectsLine(settings.landscape[i], settings.landscape[i + 1], settings.landscape[i + 2],
                        settings.landscape[i + 3])){
                    /*double[] intersection = intersection(line,
                            new Line2D.Double(new Point2D.Double(landscape[i], landscape[i + 1]),
                                    new Point2D.Double(landscape[i + 2], landscape[i + 3])));
                    distanceToLanding = getDistance(intersection);*/
                    distanceToLanding = x<settings.landingSite[0]? settings.landingSite[0]-x 
                            : x>settings.landingSite[2]? x-settings.landingSite[2] : 0;
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
            double x1 = a.getX1(), y1 = a.getY1(), x2 = a.getX2(), y2 = a.getY2(), x3 = b.getX1()
                    , y3 = b
                    .getY1(), x4 = b.getX2(), y4 = b.getY2();
            double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if (d == 0){
                return new double[]{0.0, 0.0};
            }

            double xi = ((x3 - x4) * (x1 * y2 - y1 * x2) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d;
            double yi = ((y3 - y4) * (x1 * y2 - y1 * x2) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d;

            return new double[]{xi, yi};
        }

        private boolean isSafeLand(Settings s) {
            if (x > s.landingSite[2] || x < s.landingSite[0] || angle != 0 || Math.abs(vS) > 40 ||
                Math.abs(hS) > 20) return false;
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
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            State state = (State) o;

            if (Double.compare(state.distanceToLanding, distanceToLanding) != 0) return false;
            if (Double.compare(state.x, x) != 0) return false;
            if (Double.compare(state.y, y) != 0) return false;
            if (Double.compare(state.hS, hS) != 0) return false;
            if (Double.compare(state.vS, vS) != 0) return false;
            if (fuel != state.fuel) return false;
            if (angle != state.angle) return false;
            if (power != state.power) return false;
            if (isLanded != state.isLanded) return false;
            if (isOutOfBounds != state.isOutOfBounds) return false;
            return isSafeLanded == state.isSafeLanded;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            temp   = Double.doubleToLongBits(distanceToLanding);
            result = (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(x);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(y);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(hS);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp   = Double.doubleToLongBits(vS);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            result = 31 * result + fuel;
            result = 31 * result + angle;
            result = 31 * result + power;
            result = 31 * result + (isLanded ? 1 : 0);
            result = 31 * result + (isSafeLanded ? 1 : 0);
            result = 31 * result + (isOutOfBounds ? 1 : 0);
            return result;
        }

        @Override
        public String toString() {
            return "State{" + "distanceToLanding=" + distanceToLanding + ", x=" + x + ", y=" + y +
                   ", hS=" + hS + ", vS=" + vS + ", fuel=" + fuel + ", angle=" + angle +
                   ", power=" + power + ", isLanded=" + isLanded + ", isSafeLanded=" +
                   isSafeLanded + '}';
        }
    }

    public static class Individual {
        List<Gene> genes = new ArrayList<>();
        List<State> states = new ArrayList<>();
        double fitness;

        public Individual(List<Gene> genes) {
            this.genes = genes;
        }
        public Individual() { 
            
        }
        public void removeFirstStep(){
            if(isValid()){
                genes.get(0).getValues()[2] -=1;
                if(genes.get(0).getValues()[2]<=0){
                    genes.remove(0);
                    states.remove(0);
                }
            }
        }
        public boolean isValid(){
            return !genes.isEmpty() && !states.isEmpty();
        }
        
        List<Gene> getGenes() {
            return genes;
        }

        double getFitness() {
            return fitness;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Individual that = (Individual) o;

            if (Double.compare(that.fitness, fitness) != 0) return false;
            if (genes != null ? !genes.equals(that.genes) : that.genes != null) return false;
            return states != null ? states.equals(that.states) : that.states == null;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = genes != null ? genes.hashCode() : 0;
            result = 31 * result + (states != null ? states.hashCode() : 0);
            temp   = Double.doubleToLongBits(fitness);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }

        public void simulate(State state, Settings settings) {
            int i = 0;
            states = new ArrayList<>();
            ArrayList<Gene> newGenes = new ArrayList<>();
            while(states.size()<= settings.INDIVIDUAL_LENGTH && (!state.isLanded) && (!state.isOutOfBounds)){
                Gene g;
                if(i < this.genes.size()) g = this.genes.get(i);
                else g = new Gene();
                state = state.simulate(settings,g);
                states.add(state);
                newGenes.add(g);
                i++;
            }
            genes = newGenes;
            fitness = calculateFitness( settings);
        }

        private double calculateFitness(Settings settings) {
            State lastState = states.get(states.size() - 1);
            double distance = lastState.isLanded ? lastState.distanceToLanding : 7000.0;
            double hSpeed = Math.abs(lastState.hS);
            double vSpeed = Math.abs(lastState.vS);
            int angle = Math.abs(lastState.angle);
            int fuel = Math.abs(lastState.fuel);
            //Normalize to 0-100 range by dividing current / max
            double normDistance = 10.0 - ((distance / 7000.0) * 10.0);
            //Make score exponential
            normDistance *= normDistance;
            normDistance *= normDistance;
            double normHS = 10.0 - (((hSpeed < 20.0 ? 0 : hSpeed) / 500.0) * 10.0);
            //double normHS = 10.0 - ((hSpeed / 500.0) * 10.0);
            normHS *= normHS;
            normHS *= normHS;
            double normVS = 10.0 - (((vSpeed < 40.0 ? 0 : vSpeed) / 500.0) * 10.0);
            //double normVS = 10.0 - (( vSpeed / 500.0) * 10.0);
            normVS *= normVS;
            normVS *= normVS;
            double normAngle = 10.0 - ((angle / 90.0) * 10.0);
            normAngle *= normAngle;
            normAngle *= normAngle;
            double normFuel = (fuel / (double) settings.startingFuel) * 10.0;
            normFuel *= normFuel;
            normFuel *= normFuel;

            //Get fitness
            double fitnessScore = 0;
            fitnessScore += normDistance * settings.GENE_WEIGHTS[0];
            fitnessScore += normHS * settings.GENE_WEIGHTS[1];
            fitnessScore += normVS * settings.GENE_WEIGHTS[2];
            fitnessScore += normAngle * settings.GENE_WEIGHTS[3];
            
            //Add fuel to score only after safe landing
            if (lastState.isSafeLanded) fitnessScore += normFuel * settings.GENE_WEIGHTS[4];
            return fitnessScore;
        }

        @Override
        public String toString() {
            return "Individual{" + "fitness=" +
                   fitness + "\ngenes=" + genes + "\nstates=" + states + '}';
        }
    }

    public static class Gene {
        private final int[] vars;

        public Gene() {
            vars = new int[3];
            int angleDir = Genetic.R.nextInt(3);
            vars[0] = angleDir == 0 ? -15 : angleDir == 1 ? 0 : 15;
            //vars[0] = Genetic.R.nextInt(31)-15;
            vars[1] = Genetic.R.nextInt(3) - 1;
           // vars[2] = Genetic.R.nextInt(10);
            vars[2] = 1;
        }

        public Gene(int[] vars) {
            this.vars = vars;
        }

        int[] getValues() {
            return vars;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Gene gene = (Gene) o;

            return Arrays.equals(vars, gene.vars);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(vars);
        }

        @Override
        public String toString() {
            return "Gene{" + "vars=" + Arrays.toString(vars) + '}';
        }
    }

    public static class Solver {
        public Settings settings;

        public Solver(Settings settings) {
            this.settings = settings;
        }

        public List<Individual> evolve(List<Individual> pop, State state, long startTime, int[] crossovers) {
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
                pop = newPop;
                pop.sort(Comparator.comparingDouble(Individual::getFitness).reversed());
                counter++;
            }
            crossovers[0] += counter;
            if (settings.REMOVE_STEP) removeFirstStep(pop);
            return pop;
        }
        
        private void removeFirstStep(List<Individual> pop){
            for (int i = pop.size() - 1; i >= 0; i--){
                Individual ind = pop.get(i);
                ind.removeFirstStep();
                if (!ind.isValid()) pop.remove(i);
            }
        }
        
        private void mutate(List<Individual> pop){
            for (Individual ind : pop){
                List<Gene> genes = ind.getGenes();
                for (int j = 0; j < genes.size(); j++)
                    if (Genetic.R.nextDouble() < settings.MUTATION_CHANCE) genes.set(j, new Gene());
            }
        }

        private List<Individual> crossover(List<Individual> pop) {
            List<Individual> result = new ArrayList<>();
            while (result.size() <
                   settings.DESIRED_POPULATION_SIZE * settings.CROSSOVER_PERCENTAGE) {
                Individual i1 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                Individual i2 = settings.selectType.doSelect(pop, settings.TOURNAMENT_SIZE);
                if (settings.RANDOM_CROSSOVER_ON_DUPLICATE){
                    ArrayList<Gene> randGenes = new ArrayList<>();
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
            List<Individual> result = new ArrayList<>();
            for (int i = 0; i < length; i++)
                 result.add(new Individual());
            return result;
        }
    }

    public static class Settings {
        public Genetic.SelectType selectType = SelectType.TOURNAMENT;
        public Genetic.CrossoverType crossoverType = CrossoverType.POINT;
        public boolean RANDOM_CROSSOVER_ON_DUPLICATE = true;
        public boolean REMOVE_DUPLICATES = false;
        public boolean REMOVE_STEP = false;
        public int TOURNAMENT_SIZE = 4;
        public double ELITISM_PERCENTAGE = 0.2;
        public int INDIVIDUAL_LENGTH = 200;
        public int DESIRED_POPULATION_SIZE = 50;
        public double CROSSOVER_PERCENTAGE = 0.8;
        public double MUTATION_CHANCE = 0.02;
        public double[] GENE_WEIGHTS = new double[]{0.1, 0.3, 0.4, 0.1, 1};
        public int startingFuel = 3000;
        public long EVALUATE_TIME = 90;
        public double[] landscape;
        public double[] landingSite;

        public Settings() {
        }

        public Settings(Settings s) {
            this.selectType                    = s.selectType;
            this.crossoverType                 = s.crossoverType;
            this.RANDOM_CROSSOVER_ON_DUPLICATE = s.RANDOM_CROSSOVER_ON_DUPLICATE;
            this.REMOVE_DUPLICATES             = s.REMOVE_DUPLICATES;
            this.REMOVE_STEP                   = s.REMOVE_STEP;
            this.TOURNAMENT_SIZE               = s.TOURNAMENT_SIZE;
            this.ELITISM_PERCENTAGE            = s.ELITISM_PERCENTAGE;
            this.INDIVIDUAL_LENGTH             = s.INDIVIDUAL_LENGTH;
            this.DESIRED_POPULATION_SIZE       = s.DESIRED_POPULATION_SIZE;
            this.CROSSOVER_PERCENTAGE          = s.CROSSOVER_PERCENTAGE;
            this.MUTATION_CHANCE               = s.MUTATION_CHANCE;
            this.GENE_WEIGHTS                  = s.GENE_WEIGHTS;
            this.startingFuel                  = s.startingFuel;
            this.EVALUATE_TIME                 = s.EVALUATE_TIME;
            this.landscape                     = s.landscape;
            this.landingSite                   = s.landingSite;
        }
    }
    
}



