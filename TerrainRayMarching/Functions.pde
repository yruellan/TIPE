

void setup_table(){
  //import java.io.File;
  Table table ;
  table = loadTable("../content/data.csv","header");
  for (TableRow row : table.rows()) {
    if (row.getInt("id") == 0){
      
      picture_n = row.getInt("value") ;
      println("load : picture_n =", picture_n);
    }
  }
}

void new_table(){
  Table table = new Table();
  table.addColumn("id");
  table.addColumn("name");
  table.addColumn("value");
  
  TableRow newRow = table.addRow();
  newRow.setInt("id", table.lastRowIndex());
  newRow.setString("name", "picture_n");
  newRow.setInt("value", 25);
  saveTable(table,"data2.csv");
  
}

void save_table(){
  Table table = loadTable("../content/data.csv","header");
  //print("exit() : picture_n = ",picture_n);
  for (TableRow row : table.rows()) {
    if (row.getInt("id") == 0){
      row.setInt("value",picture_n) ;
    }
  }
  saveTable(table,"../content/data.csv");
}
