function fileNames=ReadAfiFilenames(afiFileName)
    
% Convenience function to determine the svs files associated with an IF afi file
% INPUT
% afiFileName - name of the afi file
%               
% OUTPUT
% fileNames-  a cell array containing the names of the svs files corresponding to the different IF channels
    xmldata=xmlread(afiFileName);
    
    import javax.xml.xpath.*
    factory = XPathFactory.newInstance;
    xpath = factory.newXPath;
    
    % compile and evaluate the XPath Expression
    expression = xpath.compile('ImageList/Image/Path');
    nodeList= expression.evaluate(xmldata, XPathConstants.NODESET);
    
    fileNames=cell(nodeList.getLength,1);
    for i = 1:nodeList.getLength
        node = nodeList.item(i-1);
        fileNames{i}=char(node.getFirstChild.getNodeValue);
    end
end
