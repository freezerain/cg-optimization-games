package graphDrawer;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.chart.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.stage.Stage;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

public class GraphDrawerApp extends Application {
    private static final CountDownLatch latch = new CountDownLatch(1);
    public static GraphDrawerApp appReference = null;
   // private TilePane pane = new TilePane();
 //   private Stage stage;

    public GraphDrawerApp() {
        appReference = this;
    }

    public static GraphDrawerApp waitForStart() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return appReference;
    }
    
    public void showMultiLineChart(String testName, Map<String, Map<double[], Number>> data){
        TilePane pane = new TilePane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        List<LineChart<Number, Number>> chartList = new ArrayList<>();
        for (Map.Entry<String, Map<double[], Number>> var : data.entrySet()){
            String varTitle = var.getKey();
            Map<double[], Number> varMap = var.getValue();
            chartList.add(createMultiLineChart(varTitle, varMap));
        }
        pane.getChildren().addAll(chartList);
        Platform.runLater(()->{
            Stage stage = new Stage();
            stage.setWidth(1550);
            stage.setHeight(850);
            stage.setScene(new Scene(pane));
            stage.setTitle(testName);
            stage.show();
        });
    }
    private LineChart<Number, Number> createMultiLineChart(String title, Map<double[], Number> varMap){
        NumberAxis xAxis = new NumberAxis();
        NumberAxis yAxis = new NumberAxis();
        LineChart<Number, Number> lineChart = new LineChart<>(xAxis, yAxis);
        lineChart.setTitle(title);
        List<XYChart.Series<Number, Number>> seriesList = new ArrayList<>();
        for (int i = 1; i <= 5; i++){
            XYChart.Series<Number, Number> ser = new XYChart.Series<>();
            ser.setName("w"+i);
            seriesList.add(ser);
        }
        for (Map.Entry<double[], Number> e : varMap.entrySet()){
            for (int i = 0; i < 5; i++){
                seriesList.get(i).getData().add(new XYChart.Data<>(e.getKey()[i], e.getValue()));
            }
        }
        lineChart.getData().addAll(seriesList);
        return lineChart;
    }
    
    
    public void showLineChart(String testName, Map<String, Map<Number, Number>> data){
        TilePane pane = new TilePane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        List<LineChart<Number, Number>> chartList = new ArrayList<>();
        for (Map.Entry<String, Map<Number, Number>> var : data.entrySet()){
            String varTitle = var.getKey();
            Map<Number, Number> varMap = var.getValue();
            chartList.add(createLinearChart(varTitle, varMap));
        }
        pane.getChildren().addAll(chartList);
        Platform.runLater(()->{
            Stage stage = new Stage();
            stage.setWidth(1550);
            stage.setHeight(850);
            stage.setScene(new Scene(pane));
            stage.setTitle(testName);
            stage.show();
        });
    }
    
    private LineChart<Number, Number> createLinearChart(String title, Map<Number, Number> varMap){
        NumberAxis xAxis = new NumberAxis();
        NumberAxis yAxis = new NumberAxis();
        LineChart<Number, Number> lineChart = new LineChart(xAxis, yAxis);
        XYChart.Series<Number, Number> series = new XYChart.Series<>();
        series.setName(title);
        varMap.forEach((k, v) -> series.getData().add(new XYChart.Data<>(k, v)));
        lineChart.getData().add(series);
        return lineChart;
    }

    public void showBarChart(String testName,Map<String, Map<Boolean, Number>> data) {
        TilePane pane = new TilePane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        List<BarChart<String, Number>> chartList = new ArrayList<>();
        for (Map.Entry<String, Map<Boolean, Number>> var : data.entrySet()){
            String varTitle = var.getKey();
            Map<Boolean, Number> varMap = var.getValue();
            chartList.add(createBarChart(varTitle, varMap.get(true), varMap.get(false)));
        }
        pane.getChildren().addAll(chartList);
        Platform.runLater(()->{
            Stage stage = new Stage();
            stage.setWidth(1550);
            stage.setHeight(850);
            stage.setScene(new Scene(pane));
            stage.setTitle(testName);
            stage.show();
        });
    }

    private BarChart<String, Number> createBarChart(String title, Number trueVal, Number falseVal) {
        XYChart.Series<String, Number> series1 = new XYChart.Series<>();
        XYChart.Series<String, Number> series2 = new XYChart.Series<>();
        series1.setName("True");
        series2.setName("False");
        series1.getData().add(new XYChart.Data<>("", trueVal));
        series2.getData().add(new XYChart.Data<>("", falseVal));
        final BarChart<String, Number> bc = new BarChart<>(new CategoryAxis(), new NumberAxis());
        bc.setTitle(title);
        bc.getData().addAll(series1 ,series2);
        return bc;
    }

    @Override
    public void start(Stage stage) {
      //  this.stage = stage;
        /*stage.setTitle("Graph drawer");
        TilePane pane = new TilePane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        stage.setWidth(1550);
        stage.setHeight(850);
        stage.setScene(new Scene(pane));
        stage.show();*/
        latch.countDown();
    }
}
