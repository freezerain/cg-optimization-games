package marsLander;

import java.util.ArrayList;
import java.util.List;

public class LanderPath {
    List<Lander> landerList;
    Lander firstState;
    Lander lastState;
    boolean isLanded, isSafeLanded;
    double score;

    public LanderPath(List<Lander> landerList, double score) {
        this.landerList = landerList;
        firstState      = landerList.get(0);
        lastState       = landerList.get(landerList.size() - 1);
        isLanded = lastState.isLanded;
        isSafeLanded = lastState.isSafeLanded;
        this.score      = score;
    }
}
