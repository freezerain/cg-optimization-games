package src;

import javafx.application.Application;
import marsLander.MarsLanderApp;

import java.util.HashMap;
import java.util.Map;

public class DrawerController {
    int maxWidth = 1800;
    int maxHeight = 1000;
    int normalizedWidth = 1750;
    int normalizedHeight = 750;
    private MarsLanderApp app;

    public void updateCoord(int gen, double x, double y, double hSpeed, double vSpeed, int angle,
                            int power, int fuel, Map<double[], Double> paths) {
        HashMap<double[], Double> newPaths = new HashMap<double[], Double>();
        int counter = 0;
        for (Map.Entry<double[], Double> e : paths.entrySet()){
            double[] k = e.getKey();
            for (int i = 0; i < k.length; i++){
                k[i] = (i % 2 == 0 ? k[i] : 3000.0 - k[i]) / 4.0;
            }
            if(e.getValue()>=10.0) counter++;
            newPaths.put(k, e.getValue());
        }


        app.updateObservable(gen, counter,(int) Math.round(x / 4.0), (int) Math.round(750.0 - y / 4.0),
                (int) Math.round(hSpeed), (int) Math.round(vSpeed), angle, power, fuel, newPaths);
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
