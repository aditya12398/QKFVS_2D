import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;
import org.jfree.ui.*;
import java.io.*;
import javax.swing.*;

class Plot_Residual {
    private Plot_Residual(String title) throws IOException {
        plotres(title);
    }

    public void plotres(String title1) throws IOException {
        JFrame f = new JFrame(title1);
        int i, j, index;
        String line;
        Double data1[] = new Double[4];
        BufferedReader in = new BufferedReader(new FileReader("./residue"));
        XYSeriesCollection data = new XYSeriesCollection();
        XYSeries series = new XYSeries("Residue");
        while ((line = in.readLine()) != null) {
            j = 0;
            index = 0;
            for (i = 0; i < line.length(); i++) {
                if (line.charAt(i) == ' ' || line.charAt(i) == '\t' || i == (line.length() - 1)) {
                    data1[index++] = Double.parseDouble(line.substring(j, i));
                    j = i;
                }
            }
            series.add(data1[0], (Number) (Math.log(data1[2])));
        }
        data.addSeries(series);
        final JFreeChart chart = ChartFactory.createXYLineChart("Residual Plotting", "Iteration",
                "Residual Value (Log Scale)", data, PlotOrientation.VERTICAL, false, true, false);
        final ChartPanel chartPanel = new ChartPanel(chart);
        ChartUtilities.saveChartAsPNG(new File("Residue_Plot.png"), chart, 800, 570);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 570));
        f.setContentPane(chartPanel);
        f.pack();
        RefineryUtilities.centerFrameOnScreen(f);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
    }

    public static void main(String[] args) throws IOException {
        new Plot_Residual("Residue");
    }
}