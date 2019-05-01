# Installation
Step 1. Add it in your root build.gradle at the end of repositories:
```
allprojects {
		repositories {
			...
			maven { url 'https://www.jitpack.io' }
		}
	}
```  
Step 2. Add the dependency
```
dependencies {
	        implementation 'com.github.manishkummar21:AndroidPdfViewer:1.0'
	}
 ``` 
 
 # Include PDFViewActivity in your manifest file
 
 ```
 <activity android:name="com.github.pdfviewer.PDFViewActivity" />
 
 ```
 
 # Load a PDF file from anywhere in your code
 
 ```
 PDFView.with(context)
        .fromfilepath(fileapth)
        .swipeHorizontal(false) //if false pageswipe is vertical otherwise its horizontal
	.start();
```	
	
	
//Example of loading pdf file from asset folder
 ```

 File file = new File(getCacheDir(), "sample.pdf");
        if (!file.exists()) {

            try {
                InputStream asset = getAssets().open("sample.pdf");
                FileOutputStream output = null;
                output = new FileOutputStream(file);
                final byte[] buffer = new byte[1024];
                int size;
                while ((size = asset.read(buffer)) != -1) {
                    output.write(buffer, 0, size);
                }
                asset.close();
                output.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
	
	
        PDFView.with(this)
                .fromfilepath(file.getAbsolutePath()))
                .swipeHorizontal(false)
                .start()
 ```
