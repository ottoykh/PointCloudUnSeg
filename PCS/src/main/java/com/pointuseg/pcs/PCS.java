package com.pointuseg.pcs;

/**
 *
 * @author Yu Kai Him Otto
 */
import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.awt.Color;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

public class PCS {

    private TreeSet<PointCloud> points;

    public static class PointCloud implements Comparable<PointCloud> {

        float x, y, z;
        int r, g, b;

        public PointCloud(float x, float y, float z, int r, int g, int b) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.r = r;
            this.g = g;
            this.b = b;
        }

        @Override
        public int compareTo(PointCloud other) {
            return Float.compare(this.z, other.z);
        }
    }

    public class Plane {

        private final float a, b, c, d; // Coefficients of the plane equation: ax + by + cz + d = 0

        public Plane(float a, float b, float c, float d) {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
        }

        public float getA() {
            return a;
        }

        public float getB() {
            return b;
        }

        public float getC() {
            return c;
        }

        public float getD() {
            return d;
        }
    }

    private JTextField inputPathField;
    private JTextField colorField;
    private JTextField intensityValueField;
    private JTextArea messageArea;
    private JButton colorPickerButton;
    private JButton previewButton;
    private JButton segmentButton;
    private final Set<int[]> targetColors = new HashSet<>();
    private JCheckBox ransacCheckBox;
    private JCheckBox groundCheckBox;
    private Set<PointCloud> groundPoints;
    private Set<PointCloud> colorSegmentedPoints;

    public static void main(String[] args) {
        SwingUtilities.invokeLater(PCS::new);
    }

    private PCS() {
        createAndShowGUI();
    }

    private void createAndShowGUI() {
        JFrame frame = new JFrame("Point Cloud Processing");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(550, 300);

        JTabbedPane tabbedPane = new JTabbedPane();

        JPanel page1Panel = createPage1Panel(frame);
        tabbedPane.addTab("Color", page1Panel);

        JPanel page2Panel = createPage2Panel();
        tabbedPane.addTab("Intensity", page2Panel);

        frame.getContentPane().add(BorderLayout.NORTH, tabbedPane);
        frame.getContentPane().add(BorderLayout.CENTER, new JScrollPane(messageArea));

        frame.setVisible(true);
    }

    private JPanel createPage1Panel(JFrame frame) {
        JPanel panel = new JPanel(new GridLayout(4, 3));

        inputPathField = new JTextField();
        JButton browseButton = new JButton("Browse");
        colorField = new JTextField();
        colorPickerButton = new JButton("Pick Color");
        JButton processButton = new JButton("Process");
        previewButton = new JButton("Preview Point Cloud");
        messageArea = new JTextArea();
        messageArea.setEditable(false);

        groundCheckBox = new JCheckBox("Progressive Filter");
        ransacCheckBox = new JCheckBox("Perform RANSAC");

        panel.add(new JLabel(" Input path:"));
        panel.add(inputPathField);
        panel.add(browseButton);
        panel.add(new JLabel(" Pre-filter:"));
        panel.add(groundCheckBox);
        panel.add(previewButton);
        panel.add(new JLabel(" Color:"));
        panel.add(colorField);
        panel.add(colorPickerButton);
        panel.add(new JLabel(" Console Message:"));
        panel.add(ransacCheckBox);
        panel.add(processButton);

        browseButton.addActionListener(e -> browseFile());
        colorPickerButton.addActionListener(e -> pickColor());
        previewButton.addActionListener(e -> {
            try {
                String inputFilePath = inputPathField.getText();
                points = readCSVFile(inputFilePath); // Assign the loaded points to the class-level variable
                previewPointCloud(points);
            } catch (IOException ex) {
                handleIOException(ex);
                System.out.println("Error: " + ex.getMessage());
            }
        });
        processButton.addActionListener(e -> processFile());

        return panel;
    }

    private JPanel createPage2Panel() {
        JPanel panel = new JPanel(new GridLayout(4, 3));

        intensityValueField = new JTextField();
        segmentButton = new JButton("Segment");
        JButton plotHistogramButton = new JButton("Plot Histogram");
        panel.add(new JLabel(" Intensity segmentation:"));
        panel.add(plotHistogramButton);
        panel.add(new JLabel(" Intensity Threshold:"));
        panel.add(intensityValueField);
        panel.add(new JLabel(""));
        panel.add(segmentButton);
        panel.add(new JLabel(" Console Message:"));

        segmentButton.addActionListener(e -> {
            try {
                float intensityThreshold = Float.parseFloat(intensityValueField.getText());
                Map<String, TreeSet<PointCloud>> intensitySegmentedPoints = segmentPointsByIntensity(points, intensityThreshold);

                String withinThresholdFilePath = "output/Within_Threshold.pts";
                exportVerticesToCSV(intensitySegmentedPoints.get("WithinThreshold"), withinThresholdFilePath);
                previewPointCloud(intensitySegmentedPoints.get("WithinThreshold"));
                messageArea.append("[" + getCurrentTime() + "] Points within the intensity threshold exported. File: " + withinThresholdFilePath + "\n");

                // Export points outside the threshold
                String outsideThresholdFilePath = "output/Outside_Threshold.pts";
                exportVerticesToCSV(intensitySegmentedPoints.get("OutsideThreshold"), outsideThresholdFilePath);
                messageArea.append("[" + getCurrentTime() + "] Points outside the intensity threshold exported. File: " + outsideThresholdFilePath + "\n");

            } catch (NumberFormatException ex) {
                messageArea.append("Please enter a valid intensity threshold.\n");
            }
        });

        plotHistogramButton.addActionListener(e -> {
            // Use the points loaded in page 1 for the histogram
            if (points != null && !points.isEmpty()) {
                plotRHistogram(points);
            } else {
                // Handle the case where points are not loaded yet
                messageArea.append("Preview the Point Cloud first!\n");
            }
        });
        return panel;
    }

    private void pickColor() {
        Color pickedColor = JColorChooser.showDialog(null, "Choose Color", Color.BLACK);

        if (pickedColor != null) {
            targetColors.add(new int[]{pickedColor.getRed(), pickedColor.getGreen(), pickedColor.getBlue()});
            updateColorField();
        }
    }

    private void updateColorField() {
        // Convert the set of int[] to a set of Colors for display
        Set<Color> colorSet = targetColors.stream()
                .map(colorArray -> new Color(colorArray[0], colorArray[1], colorArray[2]))
                .collect(Collectors.toSet());

        StringBuilder colorStringBuilder = new StringBuilder();
        for (Color color : colorSet) {
            colorStringBuilder.append(colorToString(color)).append(", ");
        }
        colorField.setText(colorStringBuilder.toString().replaceAll(", $", ""));
    }

    private void browseFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        int result = fileChooser.showOpenDialog(null);
        if (result == JFileChooser.APPROVE_OPTION) {
            inputPathField.setText(Objects.requireNonNull(fileChooser.getSelectedFile()).getAbsolutePath());
            String inputFilePath = inputPathField.getText();

            try {
                TreeSet<PointCloud> points = readCSVFile(inputFilePath);
                findMajorColors(points);
            } catch (IOException ex) {
                handleIOException(ex);
                System.out.print("Error: " + ex.getMessage());
            }
        }
    }

    private String colorToString(Color color) {
        return String.format("%d, %d, %d", color.getRed(), color.getGreen(), color.getBlue());
    }

    private void findMajorColors(Set<PointCloud> points) {
        Map<Integer, Integer> colorCount = new HashMap<>();

        points.forEach(point -> {
            int color = point.r << 16 | point.g << 8 | point.b;
            colorCount.put(color, colorCount.getOrDefault(color, 0) + 1);
        });

        Map<Integer, Integer> sortedColors = colorCount.entrySet().stream()
                .sorted(Map.Entry.<Integer, Integer>comparingByValue().reversed())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));

        int totalPoints = points.size();
        int count = 0;

        StringBuilder majorColorsMessage = new StringBuilder("Point Cloud is loaded.");
        majorColorsMessage.append("\nMajor Colors from the Point Cloud:\n");

        for (Map.Entry<Integer, Integer> entry : sortedColors.entrySet()) {
            int color = entry.getKey();
            int r = (color >> 16) & 0xFF;
            int g = (color >> 8) & 0xFF;
            int b = color & 0xFF;

            Color textColor = new Color(r, g, b);
            String colorString = colorToString(textColor);

            majorColorsMessage.append(" Color ").append(++count).append(": ")
                    .append(colorString).append(", ")
                    .append("(").append(entry.getValue()).append("), ")
                    .append("Percentage: ").append(String.format("%.5f%%", (double) entry.getValue() / totalPoints * 100)).append("\n");

            if (count == 20) {
                break;
            }
        }
        majorColorsMessage.append("Total Points : ").append(totalPoints).append("\n");
        SwingUtilities.invokeLater(() -> {
            messageArea.setText(majorColorsMessage.toString());
        });
    }

    private void processFile() {
        SwingWorker<Void, String> worker = new SwingWorker<Void, String>() {
            @Override
            protected Void doInBackground() {
                long startTime = System.currentTimeMillis();

                String inputFilePath = inputPathField.getText();

                try {
                    publish("\n[" + getCurrentTime() + "] Process get started.");
                    TreeSet<PointCloud> points = readCSVFile(inputFilePath);
                    publish("[" + getCurrentTime() + "] Step 0: Point Cloud is imported.");

                    Files.createDirectories(Paths.get("output"));

                    // Estimate processing time for Step 1
                    estimateAndPublishProcessingTime(startTime, "Step 1");

                    float zThresholdFraction = 0.095f;
                    float xyThresholdFraction = 0.05f;

                    // Step 1: Simple/Progressive Ground Trimming
                    if (groundCheckBox.isSelected()) {
                        groundPoints = progressiveGroundTrimming(points, zThresholdFraction, xyThresholdFraction);
                        String progressiveGroundFilePath = "output/Progressive_Ground.pts";
                        exportVerticesToCSV((TreeSet<PointCloud>) groundPoints, progressiveGroundFilePath);
                        publish("[" + getCurrentTime() + "] Step 1: Progressive Ground points exported. File: " + progressiveGroundFilePath);
                    } else {
                        groundPoints = simpleGroundTrimming(points, zThresholdFraction);
                        String groundFilePath = "output/Ground.pts";
                        exportVerticesToCSV((TreeSet<PointCloud>) groundPoints, groundFilePath);
                        publish("[" + getCurrentTime() + "] Step 1: Ground points exported. File: " + groundFilePath);

                    }

                    // Step 2: Extracting Non-Ground Points
                    // Estimate processing time for Step 2
                    estimateAndPublishProcessingTime(startTime, "Step 2");
                    TreeSet<PointCloud> nonGroundPoints = extractNonGroundPoints(points, (TreeSet<PointCloud>) groundPoints);
                    String nonGroundFilePath = "output/Non_ground.pts";
                    exportVerticesToCSV(nonGroundPoints, nonGroundFilePath);
                    publish("[" + getCurrentTime() + "] Step 2: Non-ground points exported. File: " + nonGroundFilePath);

                    simulateProcessingTime();

                    // Estimate processing time for Step 3
                    estimateAndPublishProcessingTime(startTime, "Step 3");

                    // Step 3: Extracting Points by Colors
                    colorSegmentedPoints = extractPointsByColors(nonGroundPoints, targetColors);
                    String colorSegmentedFilePath = "output/Segmented.pts";
                    exportVerticesToCSV((TreeSet<PointCloud>) colorSegmentedPoints, colorSegmentedFilePath);
                    publish("[" + getCurrentTime() + "] Step 3: Color-segmented points exported. File: " + colorSegmentedFilePath);

                    // Plot top-view 3D scatter plot
                    plotTopViewScatterPlot(groundPoints, colorSegmentedPoints);

                    if (ransacCheckBox.isSelected()) {
                        // Step 4: RANSAC Segmentation with Minimum Percentage
                        // Estimate processing time for Step 4
                        estimateAndPublishProcessingTime(startTime, "Step 4");

                        double minPercentage = 0.05; // 5%
                        Set<PointCloud> ransacSegmentedPoints = ransacSegmentation((TreeSet<PointCloud>) colorSegmentedPoints, 100, 0.01, minPercentage);
                        String ransacSegmentedFilePath = "output/RANSAC_Segmented.pts";
                        exportVerticesToCSV((TreeSet<PointCloud>) ransacSegmentedPoints, ransacSegmentedFilePath);
                        publish("[" + getCurrentTime() + "] Step 4: RANSAC-segmented points exported. File: " + ransacSegmentedFilePath);

                        int ransacProgress = (int) (((float) ransacSegmentedPoints.size() / nonGroundPoints.size()) * 100);
                        publish("[" + getCurrentTime() + "] RANSAC Progress: " + ransacProgress + "%");

                    }

                } catch (IOException ex) {
                    handleIOException(ex);
                    publish("Error: " + ex.getMessage());
                }

                long endTime = System.currentTimeMillis();
                float processingTime = (float) ((endTime - startTime) * 0.001);
                publish("Total processing time: " + processingTime + " seconds.");

                return null;
            }

            @Override
            protected void process(java.util.List<String> chunks) {
                for (String message : chunks) {
                    messageArea.append(message + "\n");
                }
            }

            private void estimateAndPublishProcessingTime(long startTime, String step) {
                long currentTime = System.currentTimeMillis();
                float elapsedTime = (float) ((currentTime - startTime) * 0.001);
                float estimatedTime = estimateProcessingTime(elapsedTime);
                publish("[" + getCurrentTime() + "] Estimated processing time for " + step + ": " + estimatedTime + " seconds.");
            }
        };

        worker.execute();
    }

    private float estimateProcessingTime(float elapsedTime) {
        float constantFactor = 1.5f;
        return elapsedTime * constantFactor;
    }

    private TreeSet<PointCloud> readCSVFile(String filePath) throws IOException {
        Path path = Paths.get(filePath);
        return Files.lines(path)
                .skip(1)
                .parallel()
                .map(this::parsePointCloud)
                .filter(Objects::nonNull)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    private PointCloud parsePointCloud(String line) {
        String[] parts = line.split(" ");
        try {
            if (parts.length == 6) {
                return new PointCloud(
                        Float.parseFloat(parts[0]),
                        Float.parseFloat(parts[1]),
                        Float.parseFloat(parts[2]),
                        Integer.parseInt(parts[3]),
                        Integer.parseInt(parts[4]),
                        Integer.parseInt(parts[5])
                );
            } else if (parts.length == 4) {
                return new PointCloud(
                        Float.parseFloat(parts[0]),
                        Float.parseFloat(parts[1]),
                        Float.parseFloat(parts[2]),
                        Integer.parseInt(parts[3]), 0, 0
                );
            } else {
            }
        } catch (NumberFormatException e) {
        }
        return null;
    }

    private TreeSet<PointCloud> simpleGroundTrimming(TreeSet<PointCloud> points, float zThresholdFraction) {
        float minZ = points.first().z;
        float maxZ = points.last().z;
        float zThreshold = minZ + zThresholdFraction * (maxZ - minZ);

        return points.parallelStream().filter(point -> point.z <= zThreshold).collect(Collectors.toCollection(TreeSet::new));
    }

    private TreeSet<PointCloud> progressiveGroundTrimming(TreeSet<PointCloud> points, float zThresholdFraction, float xyThresholdFraction) {
        float minZ = points.first().z;
        float maxZ = points.last().z;
        float zThreshold = minZ + zThresholdFraction * (maxZ - minZ);

        float minX = points.stream().map(point -> point.x).min(Float::compare).orElse(0f);
        float maxX = points.stream().map(point -> point.x).max(Float::compare).orElse(0f);
        float minY = points.stream().map(point -> point.y).min(Float::compare).orElse(0f);
        float maxY = points.stream().map(point -> point.y).max(Float::compare).orElse(0f);

        float xyThresholdX = minX + xyThresholdFraction * (maxX - minX);
        float xyThresholdY = minY + xyThresholdFraction * (maxY - minY);

        return points.parallelStream()
                .filter(point -> point.z <= zThreshold && point.x >= xyThresholdX && point.y >= xyThresholdY)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    private TreeSet<PointCloud> extractNonGroundPoints(TreeSet<PointCloud> points, TreeSet<PointCloud> groundPoints) {
        TreeSet<PointCloud> nonGroundPoints = new TreeSet<>(points);
        nonGroundPoints.removeAll(groundPoints);
        return nonGroundPoints;
    }

    private TreeSet<PointCloud> extractPointsByColors(TreeSet<PointCloud> points, Set<int[]> targetColors) {
        return points.parallelStream()
                .filter(point -> targetColors.stream().anyMatch(color
                -> isWithinColorTolerance(point, color)))
                .collect(Collectors.toCollection(TreeSet::new));
    }

    private boolean isWithinColorTolerance(PointCloud point, int[] targetColor) {
        double tolerancePercentage = 0.25;

        int rDiff = Math.abs(point.r - targetColor[0]);
        int gDiff = Math.abs(point.g - targetColor[1]);
        int bDiff = Math.abs(point.b - targetColor[2]);

        int maxDiff = (int) (255 * tolerancePercentage);

        return rDiff <= maxDiff && gDiff <= maxDiff && bDiff <= maxDiff;
    }

    private void exportVerticesToCSV(TreeSet<PointCloud> points, String filePath) {
        Path path = Paths.get(filePath);
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(path))) {
            points.parallelStream().forEach(point
                    -> writer.format("%.4f %.4f %.4f %d %d %d%n", point.x, point.y, point.z, point.r, point.g, point.b));
        } catch (IOException e) {
            handleIOException(e);
        }
    }

    private void handleIOException(IOException e) {
        e.printStackTrace();
    }

    private String getCurrentTime() {
        return LocalDateTime.now().format(DateTimeFormatter.ofPattern("HH:mm:ss"));
    }

    private void simulateProcessingTime() {
        try {
            TimeUnit.SECONDS.sleep(2);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            handleIOException(new IOException("Thread sleep interrupted", e));
        }
    }

    private Set<PointCloud> ransacSegmentation(TreeSet<PointCloud> colorSegmentedPoints, int iterations, double distanceThreshold, double minPercentage) {
        Set<PointCloud> segmentedPoints = new HashSet<>();
        int minPoints = (int) (colorSegmentedPoints.size() * minPercentage);

        PointCloud[] randomPoints = new PointCloud[3];

        for (int i = 0; i < iterations; i++) {
            Iterator<PointCloud> iterator = colorSegmentedPoints.iterator();
            for (int j = 0; j < 3; j++) {
                randomPoints[j] = iterator.next();
            }

            Plane fittedPlane = fitPlane(randomPoints[0], randomPoints[1], randomPoints[2]);

            Set<PointCloud> currentInliers = new HashSet<>();
            for (PointCloud point : colorSegmentedPoints) {
                if (calculateDistanceToPlane(point, fittedPlane) < distanceThreshold) {
                    currentInliers.add(point);
                }
            }

            if (currentInliers.size() > segmentedPoints.size() && currentInliers.size() >= minPoints) {
                segmentedPoints.clear();
                segmentedPoints.addAll(currentInliers);
            }
        }

        return segmentedPoints;
    }

    private Plane fitPlane(PointCloud p1, PointCloud p2, PointCloud p3) {
        float x1 = p1.x, y1 = p1.y, z1 = p1.z;
        float x2 = p2.x, y2 = p2.y, z2 = p2.z;
        float x3 = p3.x, y3 = p3.y, z3 = p3.z;

        float a = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
        float b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
        float c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
        float d = -a * x1 - b * y1 - c * z1;

        return new Plane(a, b, c, d);
    }

    private double calculateDistanceToPlane(PointCloud point, Plane plane) {
        float x = point.x;
        float y = point.y;
        float z = point.z;

        float a = plane.getA();
        float b = plane.getB();
        float c = plane.getC();
        float d = plane.getD();

        return Math.abs(a * x + b * y + c * z + d) / Math.sqrt(a * a + b * b + c * c);
    }

    private void previewPointCloud(Set<PointCloud> points) {
        SwingUtilities.invokeLater(() -> {
            JFrame previewFrame = new JFrame("Preview Point Cloud");
            previewFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            previewFrame.setSize(600, 600);

            JPanel previewPanel = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Graphics2D g2d = (Graphics2D) g;
                    renderTopViewScatterPlot(g2d, getWidth(), getHeight(), points, this);
                }
            };

            previewFrame.getContentPane().add(previewPanel);
            previewFrame.setVisible(true);
        });
    }

    private void renderTopViewScatterPlot(Graphics2D g2d, int width, int height, Set<PointCloud> points, JPanel panel) {
        if (points != null && !points.isEmpty()) {
            float minX = Float.MAX_VALUE, minY = Float.MAX_VALUE, maxX = Float.MIN_VALUE, maxY = Float.MIN_VALUE;

            for (PointCloud point : points) {
                minX = Math.min(minX, point.x);
                minY = Math.min(minY, point.y);
                maxX = Math.max(maxX, point.x);
                maxY = Math.max(maxY, point.y);
            }

            float xScale = width / (maxX - minX);
            float yScale = height / (maxY - minY);

            for (PointCloud point : points) {
                int x = (int) ((point.x - minX) * xScale);
                int y = (int) ((point.y - minY) * yScale);
                int size = 5;
                g2d.setColor(new Color(point.r, point.g, point.b)); // Set color to RGB values
                g2d.fillOval(x, y, size, size);
            }

            panel.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    float clickedX = e.getX();
                    float clickedY = e.getY();
                    float minX = Float.MAX_VALUE, minY = Float.MAX_VALUE, maxX = Float.MIN_VALUE, maxY = Float.MIN_VALUE;

                    for (PointCloud point : points) {
                        minX = Math.min(minX, point.x);
                        minY = Math.min(minY, point.y);
                        maxX = Math.max(maxX, point.x);
                        maxY = Math.max(maxY, point.y);
                    }

                    float xScale = width / (maxX - minX);
                    float yScale = height / (maxY - minY);

                    float clickedXValue = (clickedX / xScale) + minX;
                    float clickedYValue = (clickedY / yScale) + minY;
                    float clickedZValue = findClickedZ(points, clickedXValue, clickedYValue);

                    String message = "[" + getCurrentTime() + "] Point: X=" + clickedXValue + ", Y=" + clickedYValue + ", Z=" + clickedZValue;

                    SwingUtilities.invokeLater(() -> {
                        messageArea.append(message + "\n");
                    });
                }
            });
        }
    }

    private void plotRHistogram(Set<PointCloud> points) {
        SwingUtilities.invokeLater(() -> {
            JFrame histogramFrame = new JFrame("Intensity Histogram");
            histogramFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            histogramFrame.setSize(600, 280);

            JPanel histogramPanel = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Graphics2D g2d = (Graphics2D) g;
                    renderRHistogram(g2d, getWidth(), getHeight(), (TreeSet<PointCloud>) points, this);
                }
            };

            histogramFrame.getContentPane().add(histogramPanel);
            histogramFrame.setVisible(true);
        });
    }

    private void renderRHistogram(Graphics2D g2d, int width, int height, TreeSet<PointCloud> points, JPanel panel) {
        if (points != null && !points.isEmpty()) {
            float maxIntensity = points.stream().map(point -> point.r).max(Float::compare).orElse(0);
            int[] intensityCounts = new int[256];

            for (PointCloud point : points) {
                int intensityBin = Math.round((point.r / maxIntensity) * 255);
                intensityCounts[intensityBin]++;
            }

            int maxCount = Arrays.stream(intensityCounts).max().orElse(1);

            int barWidth = width / 256;
            g2d.setColor(Color.RED);

            for (int i = 0; i < intensityCounts.length; i++) {
                int barHeight = (int) (((double) intensityCounts[i] / maxCount) * height);
                int x = i * barWidth;
                int y = height - barHeight;
                g2d.fillRect(x, y, barWidth, barHeight);
            }

            g2d.setColor(Color.BLACK);
            g2d.drawLine(0, height, width, height); // X-axis
            g2d.drawLine(0, 0, 0, height);         // Y-axis

            panel.addMouseListener(new MouseAdapter() {
                @Override
                public void mouseClicked(MouseEvent e) {
                    int clickedX = e.getX();
                    int binWidth = width / 256;
                    int clickedBin = clickedX / binWidth;

                    if (clickedBin >= 0 && clickedBin < intensityCounts.length) {
                        int intensityValue = Math.round((float) clickedBin / 255 * maxIntensity);
                        String message = "Intensity value: " + intensityValue + ", Count :" + clickedBin;
                        SwingUtilities.invokeLater(() -> {
                            messageArea.append(message + "\n");
                        });
                    }
                }
            });

            for (int i = 0; i <= 255; i += 64) {
                int x = i * (width / 255);
                g2d.drawString(Integer.toString(i), x, height + 15);
            }
        }
    }

    private Map<String, TreeSet<PointCloud>> segmentPointsByIntensity(TreeSet<PointCloud> points, float intensityThreshold) {
        Map<String, TreeSet<PointCloud>> segmentedPoints = new HashMap<>();
        TreeSet<PointCloud> withinThreshold = new TreeSet<>();
        TreeSet<PointCloud> outsideThreshold = new TreeSet<>();

        for (PointCloud point : points) {
            if (point.r >= intensityThreshold) {
                withinThreshold.add(point);
            } else {
                outsideThreshold.add(point);
            }
        }

        segmentedPoints.put("WithinThreshold", withinThreshold);
        segmentedPoints.put("OutsideThreshold", outsideThreshold);

        return segmentedPoints;
    }

    private float findClickedZ(Set<PointCloud> points, float clickedX, float clickedY) {
        PointCloud closestPoint = null;
        float minDistance = Float.MAX_VALUE;

        for (PointCloud point : points) {
            float distance = (float) Math.sqrt(Math.pow(point.x - clickedX, 2) + Math.pow(point.y - clickedY, 2));
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = point;
            }
        }

        if (closestPoint != null) {
            return closestPoint.z;
        }

        return 0.0f;
    }

    private void plotTopViewScatterPlot(Set<PointCloud> groundPoints, Set<PointCloud> segmentedPoints) {
        SwingUtilities.invokeLater(() -> {
            JFrame groundFrame = new JFrame("Ground Points (Top View)");
            groundFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            groundFrame.setSize(600, 600);

            JFrame segmentedFrame = new JFrame("Segmented Points (Top View)");
            segmentedFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            segmentedFrame.setSize(600, 600);

            JPanel groundPanel = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Graphics2D g2d = (Graphics2D) g;
                    renderTopViewScatterPlot(g2d, getWidth(), getHeight(), groundPoints, this);
                }
            };

            JPanel segmentedPanel = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Graphics2D g2d = (Graphics2D) g;
                    renderTopViewScatterPlot(g2d, getWidth(), getHeight(), segmentedPoints, this);
                }
            };

            groundFrame.getContentPane().add(groundPanel);
            segmentedFrame.getContentPane().add(segmentedPanel);

            groundFrame.setVisible(true);
            segmentedFrame.setVisible(true);
        });
    }
}
