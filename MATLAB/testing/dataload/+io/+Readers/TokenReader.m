classdef TokenReader < io.Readers.AbstractReader

    properties( Dependent = true )

        CurrentToken 
        TokenDescription    (1,:) char;
        NumberValue         (1,1) double
        WordValue           (1,:) char
        
        StoredNumbers       (:,1) cell
        StoredWords         (:,1) cell
        StoredCharacters    (:,1) cell;
        StoredComments      (:,1) cell;
        
        LineNumber          (1,1) double
        
        StoreResults        (1,1) logical
        StoreByLine         (1,1) logical
        
        EoLIsSignificant    (1,1) logical
        
        ParseNumbers        (1,1) logical
        
        LowerCaseMode       (1,1) logical
    end

    properties(Access = protected)

        NUMBERS_            (:,1) cell
        WORDS_              (:,1) cell
        CHARACTERS_         (:,1) cell;
        COMMENTS_           (:,1) cell;
        
        currentLine         (1,1) uint32 = 1;
        
        lineNumberBuffer    (1,:) cell
        lineWordBuffer      (1,:) cell
        lineCharBuffer      (1,:) cell
        lineCommentBuffer   (1,:) cell        
        
        StoreResults_       (1,1) logical = false;
        StoreByLine_        (1,1) logical = true;
        
        EoLIsSignificant_   (1,1) logical = true;  
        
        ParseNumbers_       (1,1) logical = true
        
        LowerCaseMode_      (1,1) logical = false;
    end

    properties(Access = public, Constant = true)

        EOF      = java.io.StreamTokenizer.TT_EOF;
        EOL      = java.io.StreamTokenizer.TT_EOL;
        NUMBER   = java.io.StreamTokenizer.TT_NUMBER;
        WORD     = java.io.StreamTokenizer.TT_WORD;
        INIT     = int8(-4);
    end
    
    methods
        function this = TokenReader( entry, inputStream, buffSize )

            this@io.Readers.AbstractReader( entry, inputStream, buffSize );

        end
    end 
    
    methods
        
        function numbers = get.StoredNumbers( this )
           
            if this.StoreResults_
                numbers = this.NUMBERS_;
            else
                numbers = [];
            end
        end

        function words = get.StoredWords( this )
           
            if this.StoreResults_
                words = this.WORDS_;
            else
                words = {};
            end
        end 
        
        function characters = get.StoredCharacters( this )
           
            if this.StoreResults_
                characters = this.CHARACTERS_;
            else
                characters = {};
            end
        end
        
        function comments = get.StoredComments( this )
           
            if this.StoreResults_
                comments = this.COMMENTS_;
            else
                comments = {};
            end
        end   
        
        function token = get.CurrentToken( this )
            
            ttype = this.Reader.ttype;
            if (ttype < 0) || isequal(ttype,10)
                token = ttype;
            else
                token = char(ttype);
            end
        end

        function desc = get.TokenDescription( this )
            ttype = this.Reader.ttype;
            switch ttype
                case this.EOF
                    desc = 'End of file';
                 case this.EOL
                    desc = 'End of line';                   
                 case this.NUMBER
                    desc = 'Number';
                 case this.WORD
                    desc = 'Word';
                 case this.INIT
                    desc = 'Unread';  
                otherwise
                    desc = 'Single character or quoted text';
                
            end
        end
        
        function nval = get.NumberValue( this )
            
            nval = this.Reader.nval;
        end
        
        function sval = get.WordValue( this )
            
            sval = char(this.Reader.sval);
            
            ttype = this.Reader.ttype;
            if (ttype > 0) && (ttype ~= 10)
                sval = char(ttype);
            end
        end 
        
        function lineno = get.LineNumber( this )
           
            lineno = this.Reader.lineno;
        end
        
        function val = get.StoreResults( this )
            
           val =  this.StoreResults_;
        end
        
        function set.StoreResults( this, val )
            
            this.StoreResults_ = val;
        end
        
        function val = get.StoreByLine( this )
            
           val =  this.StoreByLine_;
        end 
        
        function set.StoreByLine( this, val )
            
            this.StoreByLine_ = val;
        end        
        
        function val = get.EoLIsSignificant( this )
            
           val =  this.EoLIsSignificant_;
        end
        
        function set.EoLIsSignificant( this, val )
            
            this.validateBoolean(val,'set.EoLIsSignificant')
            
            this.Reader.eolIsSignificant(val);            
            
            this.EoLIsSignificant_ = val;
        end
        
        function val = get.LowerCaseMode( this )
            
            val = this.LowerCaseMode_;
        end
        
        function set.LowerCaseMode( this, val )
            
            this.validateBoolean(val,'set.LowerCaseMode')
            
            this.Reader.lowerCaseMode(val);
            
            this.LowerCaseMode_ = val;
        end
        
        function val = get.ParseNumbers( this )
            
           val =  this.ParseNumbers_;
        end
        
        function set.ParseNumbers( this, val )
            
            %TODO: set java reader
            this.ParseNumbers_ = logical(val);
        end        
    end
    
    methods(Access = public)
        
        function this = parseAll( this )
            
            this.StoreResults = true;
            
            isEoF = false;
            
            while ~isEoF
               isEoF = this.nextToken(); 
            end
            
        end
        
        function isEoF = nextToken( this )
            
            try
                this.Reader.nextToken();
            catch ME
               io.Util.handleIOException(ME); 
            end
            
            ttype = this.Reader.ttype;
            
            isEoF = (ttype == this.EOF);
            
            if isEoF
                return
            end
            
            if this.StoreResults
                
                if this.StoreByLine
                    
                    switch ttype
                        case this.EOL
                            ln = this.Reader.lineno - 1;
                            this.NUMBERS_{ln,:} = this.lineNumberBuffer;
                            this.WORDS_{ln,:} = this.lineWordBuffer;
                            this.CHARACTERS_{ln,:} = this.lineCharBuffer;
                            this.COMMENTS_{ln,:} = this.lineCommentBuffer;                            
                            this.lineNumberBuffer = {};
                            this.lineWordBuffer = {};
                            this.lineCharBuffer = {};
                            this.lineCommentBuffer = {};
                        case this.NUMBER
                            number = double(this.Reader.nval);
                            try this.lineNumberBuffer{1,end+1} = number; catch; end
                        case this.WORD
                            word = char(this.Reader.sval);
                            try this.lineWordBuffer{1,end+1} = word; catch; end
                        otherwise
                            word = char(ttype);
                            if isequal(length(word),1)
                                try this.lineCharBuffer{1,end+1} = word; catch; end
                            else
                                try this.lineCommentBuffer{1,end+1} = word; catch; end
                            end
                    end

                else
                
                    switch ttype
                        case this.NUMBER
                            number = double(this.Reader.nval);
                            if isempty(this.NUMBERS_)
                                this.NUMBERS_{1} = number;
                            else
                                this.NUMBERS_{end+1} = number;
                            end
                        case this.WORD
                            word = char(this.Reader.sval);
                            if isempty(this.WORDS_)
                                this.WORDS_{1} = word;
                            else
                                this.WORDS_{end+1} = word;
                            end
                    end
                end
            end
        end

        function commentChar( this, ch )
            this.Reader.commentChar(int8(ch));
        end
        
        function lineno = lineno( this )
           lineno = this.Reader.lineno(); 
        end

        function ordinaryChar( this, ch )
            this.Reader.ordinaryChar(int8(ch));
        end

        function ordinaryChars( this, chlow, chhigh )
            this.Reader.ordinaryChars(int8(chlow),int8(chhigh));
        end  

        function parseNumbers( this )
           this.Reader.parseNumbers(); 
        end

        function pushBack( this )
            this.Reader.pushBack(); 
        end

        function quoteChar( this, ch )
            this.Reader.quoteChar(int8(ch));
        end

        function resetSyntax( this )
            this.Reader.resetSyntax(); 
        end 

        function slashSlashComments( this, flag )
            this.Reader.slashSlashComments(logical(flag));
        end 

        function slashStarComments( this, flag )
            this.Reader.slashStarComments(logical(flag));
        end 

        function whitespaceChars( this, chlow, chhigh )
            this.Reader.whitespaceChars(int8(chlow),int8(chhigh));
        end 

        function wordChars( this, chlow, chhigh )
            this.Reader.wordChars(int8(chlow),int8(chhigh));
        end
        
    end

    methods(Access = protected)

        function createReader( this, inputStream, ~ )

            this.Reader = java.io.StreamTokenizer(...
                java.io.BufferedReader(...
                java.io.InputStreamReader(inputStream)));
            
            this.Reader.eolIsSignificant(true);
            
            this.Reader.parseNumbers();
        end

    end
end

