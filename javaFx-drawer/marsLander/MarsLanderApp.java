package marsLander;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Polyline;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.stage.Stage;

import java.util.*;
import java.util.concurrent.CountDownLatch;

public class MarsLanderApp extends Application {
    private static final CountDownLatch latch = new CountDownLatch(1);
    public static MarsLanderApp appReference = null;
   
    int WIDTH = 1750;
    int HEIGHT = 750;
    int HALF_IMAGE_WIDTH = 33;
    int HALF_IMAGE_HEIGHT = 31;
    private Pane pane;
    private List<Polyline> previousPaths = new ArrayList<>();
    private final SimpleIntegerProperty y = new SimpleIntegerProperty();
    private final SimpleIntegerProperty x = new SimpleIntegerProperty();
    private final SimpleIntegerProperty rot = new SimpleIntegerProperty();
    private final SimpleStringProperty coordText = new SimpleStringProperty();
    private final SimpleStringProperty speedText = new SimpleStringProperty();
    private final SimpleStringProperty controlText = new SimpleStringProperty();
    

    public MarsLanderApp() {
        appReference = this;
        latch.countDown();
    }

    public static void main(String[] args) {
        launch();
    }

    public static MarsLanderApp waitForStart() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return appReference;
    }

    public void launchApp() {
        launch();
    }

    @Override
    public void start(Stage stage) {
        stage.setTitle("Mars Lander");
        pane = new Pane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        Canvas canvas = new CanvasGridGenerator(WIDTH, HEIGHT, 10).generate(0.25,
                Color.color(1, 1, 1, 0.01));
        pane.getChildren().add(canvas);
        
        Text text1 = new Text();
        text1.setFill(Color.WHITE);
        text1.setFont(Font.font(18.0));
        text1.textProperty().bind(coordText);
        Text text2 = new Text();
        text2.setFill(Color.WHITE);
        text2.setFont(Font.font(18.0));
        text2.textProperty().bind(speedText);
        Text text3 = new Text();
        text3.setFill(Color.WHITE);
        text3.setFont(Font.font(18.0));
        text3.textProperty().bind(controlText);
        VBox vbox = new VBox(text1, text2,text3);
        vbox.setLayoutX(100);
        vbox.setLayoutY(100);
        pane.getChildren().add(vbox);
        
        Image image = new Image(getClass().getResource("MarsLander1.png").toString());
        ImageView im = new ImageView(image);
        im.layoutXProperty().bind(y);
        im.layoutYProperty().bind(x);
        im.rotateProperty().bind(rot);
        pane.getChildren().add(im);
        stage.setWidth(WIDTH);
        stage.setHeight(HEIGHT + 20);
        stage.setScene(new Scene(pane));
        stage.show();
    }

    public void updateObservable(int x, int y, int hSpeed, int vSpeed, int angle, int power, int fuel, Map<double[], Double> paths) {
        this.y.set(x - HALF_IMAGE_WIDTH);
        this.x.set(y - HALF_IMAGE_HEIGHT);
        rot.set(-angle);
        coordText.set("coord x: " + x + " y: " + y);
        speedText.set("speed x: " + hSpeed + " y: " + vSpeed);
        controlText.set("angle: " + angle + " power: " + power + " fuel: " + fuel);
        buildPaths(paths);
    }

    private void buildPaths(Map<double[], Double> paths) {
        Double maxScore = paths.entrySet()
                .stream()
                .max(Comparator.comparingDouble(Map.Entry::getValue))
                .map(Map.Entry::getValue).orElse(1.0);
        List<Polyline> newPath = new ArrayList<>();
        for (Map.Entry<double[], Double> e : paths.entrySet()){
            double[] path = e.getKey();
            Polyline polyline = new Polyline(path);
            polyline.setStroke(Color.color(1.0-(e.getValue()/maxScore), e.getValue()/maxScore, 0));
            newPath.add(polyline);
        }
        Platform.runLater(() -> {
            pane.getChildren().removeAll(previousPaths);
            pane.getChildren().addAll(newPath);
            previousPaths = newPath;
        });
    }

    public void buildLandscape(double[] landscape) {
        Polyline land = new Polyline(landscape);
        land.setStrokeWidth(1);
        land.setStroke(Color.DARKRED);
        pane.getChildren().add(land);
    }
}
