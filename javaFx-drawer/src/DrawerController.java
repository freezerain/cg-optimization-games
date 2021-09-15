package src;

import javafx.application.Application;
import marsLander.LanderPath;
import marsLander.MarsLanderApp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DrawerController {
    int maxWidth = 1800;
    int maxHeight = 1000;
    int normalizedWidth = 1750;
    int normalizedHeight = 750;
    private MarsLanderApp app;

    public void updateCoord(int gen, List<LanderPath> paths) {
        app.updateObservable(gen, paths);
    }

    public void startApp(double[] landscape) {
        long time = System.currentTimeMillis();
        new Thread(() -> Application.launch(MarsLanderApp.class)).start();
        double[] normLand = new double[landscape.length];
        for (int i = 0; i < landscape.length; i++){
            if (i % 2 != 0) normLand[i] = 750.0 - landscape[i] / 4.0;
            else normLand[i] = landscape[i] / 4.0;
        }
        app = MarsLanderApp.waitForStart();
        app.buildLandscape(normLand);
        System.out.println("await finished: " + (System.currentTimeMillis() - time) + "ms.");
    }

}
